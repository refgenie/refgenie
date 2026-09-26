/**
 * The build form.
 *
 * A full page rather than a modal: the field count is recipe-dependent and can
 * exceed a modal's usable height, and modals are reserved for confirmations and
 * short forms.
 *
 * Form state is `useState` plus a small reducer for the dynamic field map. No
 * form library: the schema is dynamic and per-recipe, which is the case those
 * libraries handle worst, and a dependency added for one form is a dependency
 * to maintain forever.
 */

import { useCallback, useEffect, useMemo, useRef, useState } from 'react';
import { useNavigate, useSearchParams } from 'react-router-dom';
import { useLocalApiClient } from '../hooks/useApiClient';
import { useCapability } from '../hooks/useCapability';
import { useConfigurations } from '../hooks/queries/useConfigurations';
import { useToast } from '../hooks/useToast';
import { buildAsset, preflightBuild } from '../services/actions';
import { ApiError } from '../services/http';
import { provisionalJob, useJobStore } from '../stores/jobStore';
import { Breadcrumbs } from '../components/common/Breadcrumbs';
import { MiniHero } from '../components/layout/MiniHero';
import { FormField } from '../components/common/FormField';
import { useHelpDisclosure } from '../components/common/useHelpDisclosure';
import { GenomeSelect } from '../components/build/GenomeSelect';
import { RecipeSelect } from '../components/build/RecipeSelect';
import { useRecipeCatalog } from '../components/build/recipeCatalog';
import { RecipeInputFields } from '../components/build/RecipeInputFields';
import { InitGenomeModal } from '../components/actions/InitGenomeModal';
import { NotAvailablePage } from './NotAvailablePage';
import {
  emptyBuildForm,
  registryPathPreview,
  validateBuildForm,
} from '../components/build/validate';
import type { BuildFormValues, FieldErrors } from '../components/build/validate';
import type { InputSection } from '../components/build/RecipeInputFields';
import type { BuildRequest } from '../services/contracts';
import type { RecipePublic } from '../types/api';

const PREFLIGHT_DEBOUNCE_MS = 400;

/**
 * Form state -> the wire body.
 *
 * Two shape facts the backend enforces with `extra="forbid"`: the recipe
 * inputs are NESTED under `params`, and there is no `docker` field at all.
 */
function toRequest(values: BuildFormValues): BuildRequest {
  return {
    recipe: values.recipeName ?? '',
    genome: values.genome,
    asset_group: values.assetGroupName,
    asset: values.assetName || null,
    recipe_version: values.recipeVersion,
    description: values.assetDescription || null,
    stage: values.stage,
    pull_parents: values.pullParents,
    params: {
      params: values.params,
      files: values.files,
      assets: values.assets,
    },
  };
}

export function BuildPage() {
  const canBuild = useCapability('build');
  const canInitGenome = useCapability('genome_init');
  const client = useLocalApiClient();
  const toast = useToast();
  const navigate = useNavigate();
  const [params] = useSearchParams();
  const catalog = useRecipeCatalog();

  const [values, setValues] = useState<BuildFormValues>(emptyBuildForm);
  const [serverErrors, setServerErrors] = useState<FieldErrors>({});
  /*
   * Preflight's `asset_name_required`: the recipe names its output after a
   * tool version, and reading that version came back empty. That is a field
   * refgenie tried to fill in and could not, not a failure -- the same request
   * with a name preflights ok -- so it is held apart from `serverErrors` and
   * turns Asset name into a required field instead of a red block of prose.
   * Holds the server's explanation, which goes behind the field's "?".
   */
  const [nameRequest, setNameRequest] = useState<string | null>(null);
  const [resolvedPath, setResolvedPath] = useState<string | null>(null);
  const [submitting, setSubmitting] = useState(false);
  const [attempted, setAttempted] = useState(false);
  const [initOpen, setInitOpen] = useState(false);
  const prefilled = useRef(false);
  /*
   * Whether the asset group in the box is the user's or the recipe's. An
   * untouched group follows the recipe, so switching recipes does not leave
   * `epilog_index` sitting over a bowtie2 build; a typed one survives, because
   * silently overwriting what someone typed is worse than a stale default.
   * Emptying the box hands the field back to the recipe.
   */
  const groupEdited = useRef(false);
  // The staging checkbox gets the same "?" affordance as a FormField, which it
  // cannot inherit: a checkbox has no FormField wrapper to collapse its hint.
  const stageHelp = useHelpDisclosure('build-stage-hint', 'staging');

  // Staging without a stage folder raises a bare ValueError server-side, so the
  // checkbox is only offered when the configuration actually has one.
  const configurations = useConfigurations({ limit: 1 });
  const stageFolder = configurations.data?.items?.[0]?.genome_stage_folder ?? null;

  const recipe: RecipePublic | null = useMemo(() => {
    if (!values.recipeName) return null;
    const versions = catalog.byName.get(values.recipeName) ?? [];
    return versions.find((r) => r.version === values.recipeVersion) ?? versions[0] ?? null;
  }, [catalog, values.recipeName, values.recipeVersion]);

  // Prefill once, from ?genome=&recipe=&asset_group=, so "build this" links work.
  useEffect(() => {
    if (prefilled.current || catalog.isPending) return;
    prefilled.current = true;
    // `genome` is an ALIAS: preflight resolves it with `alias.resolve`, which
    // does not accept a digest. The digest rides along separately, for the job
    // target key only.
    const genome = params.get('genome');
    const genomeDigest = params.get('genome_digest');
    const recipeName = params.get('recipe');
    const assetGroup = params.get('asset_group');
    if (!genome && !genomeDigest && !recipeName && !assetGroup) return;
    // An explicit `?asset_group=` is the caller's choice, not a default.
    if (assetGroup) groupEdited.current = true;
    const versions = recipeName ? (catalog.byName.get(recipeName) ?? []) : [];
    setValues((current) => ({
      ...current,
      genome: genome ?? current.genome,
      genomeDigest: genomeDigest ?? current.genomeDigest,
      recipeName: versions[0]?.name ?? current.recipeName,
      recipeVersion: versions[0]?.version ?? current.recipeVersion,
      assetGroupName: assetGroup ?? versions[0]?.name ?? current.assetGroupName,
    }));
  }, [catalog, params]);

  const clientErrors = useMemo(
    () => validateBuildForm(recipe, values, { assetNameRequired: nameRequest !== null }),
    [recipe, values, nameRequest],
  );
  const errors: FieldErrors = { ...serverErrors, ...clientErrors };
  const valid = Object.keys(clientErrors).length === 0;
  // The blank-name-while-required error blocks the SUBMIT, not the preflight:
  // preflight is what clears the requirement once a name is typed, and it
  // must keep checking the other inputs while the name box is still empty.
  const preflightable = useMemo(
    () => Object.keys(validateBuildForm(recipe, values)).length === 0,
    [recipe, values],
  );

  const setField = useCallback(<K extends keyof BuildFormValues>(key: K, value: BuildFormValues[K]) => {
    setValues((current) => ({ ...current, [key]: value }));
  }, []);

  const setInput = useCallback((section: InputSection, name: string, value: string) => {
    setValues((current) => ({ ...current, [section]: { ...current[section], [name]: value } }));
  }, []);

  // Preflight is the only thing that can tell the user a path does not exist or
  // a parent asset is missing, so it runs on a debounce once the form is
  // otherwise valid — never as a blocking round-trip on every keystroke.
  useEffect(() => {
    if (!preflightable || !canBuild) return;
    const timer = setTimeout(() => {
      preflightBuild(client, toRequest(values))
        .then((response) => {
          const next: FieldErrors = {};
          let askForName: string | null = null;
          for (const item of response.errors ?? []) {
            if (item.field === 'asset' && item.code === 'asset_name_required') {
              askForName = item.message;
            } else if (item.field) {
              next[item.field] = item.message;
            }
          }
          setServerErrors(next);
          setNameRequest(askForName);
          // `resolved` has no registry_path: compose it from what it does say.
          const assetName = response.resolved?.asset_name;
          setResolvedPath(
            assetName ? `${values.genome}/${values.assetGroupName}:${assetName}` : null,
          );
        })
        .catch(() => {
          // Preflight is advisory. If it is unavailable the submit still works
          // and the job reports the failure instead.
          setServerErrors({});
          setNameRequest(null);
        });
    }, PREFLIGHT_DEBOUNCE_MS);
    return () => clearTimeout(timer);
  }, [client, values, preflightable, canBuild]);

  if (!canBuild) {
    return (
      <NotAvailablePage
        title="Build"
        reason="This instance does not expose building (capability `build` is off)."
      />
    );
  }

  const handleRecipeChange = (next: RecipePublic | null) => {
    setValues((current) => ({
      ...current,
      recipeName: next?.name ?? null,
      recipeVersion: next?.version ?? null,
      // The CLI's `refgenie build genome/asset_group` convention: the group
      // defaults to the recipe name.
      assetGroupName: groupEdited.current ? current.assetGroupName : (next?.name ?? ''),
      params: {},
      files: {},
      assets: {},
    }));
    setServerErrors({});
    setNameRequest(null);
  };

  const handleSubmit = async (event: React.FormEvent) => {
    event.preventDefault();
    setAttempted(true);
    if (!valid) return;
    setSubmitting(true);
    try {
      const ref = await buildAsset(client, toRequest(values));
      const store = useJobStore.getState();
      if (ref.duplicate) {
        store.focusJob(ref.job_id);
        toast.info('That build is already running.');
      } else {
        store.registerQueued(
          provisionalJob({
            id: ref.job_id,
            kind: ref.kind ?? 'build',
            status: ref.status ?? 'queued',
            created_at: ref.created_at,
            label: `build ${registryPathPreview(values)}`,
            target: {
              genome_digest: values.genomeDigest,
              genome_name: values.genome,
              asset_group_name: values.assetGroupName,
              asset_name: values.assetName || null,
            },
          }),
        );
        toast.success('Build queued. Progress is in the job console.');
      }
      // The build's feedback lives in the console, not on this form.
      navigate(values.genomeDigest ? `/genomes/${values.genomeDigest}` : '/genomes');
    } catch (error) {
      if (error instanceof ApiError && error.field) {
        setServerErrors((current) => ({ ...current, [error.field as string]: error.detail }));
      }
      toast.error(error instanceof ApiError ? error.detail : 'Could not submit the build.');
    } finally {
      setSubmitting(false);
    }
  };

  const show = (key: string) => (attempted ? errors[key] : serverErrors[key]);

  /*
   * Preflight also reports coarse fields that no single input owns — `stage`,
   * `params`, and `params.assets` when the whole resolution failed. Rendering
   * only the fine-grained ones would silently drop a real reason the build
   * cannot start, so anything unclaimed surfaces above the submit button.
   */
  const claimsAnInput = (field: string) =>
    ['genome', 'recipe', 'asset_group', 'asset'].includes(field) ||
    /^params\.(params|files|assets)\./.test(field);
  const formLevelErrors = Object.entries(serverErrors).filter(
    ([field]) => !claimsAnInput(field),
  );

  return (
    // The measure cap belongs on the FORM, never around the MiniHero. The head
    // full-bleeds itself with a negative inline margin of `50% - 50vw`, and
    // that 50% is half of its PARENT -- so a parent narrower than the content
    // column drags the head off the left edge of the screen. This page's
    // 768px cap put the <h1> at x = -164 on a 1400px viewport.
    <div className="flex flex-col gap-6">
      <MiniHero
        title="Build asset"
        documentTitle="Build"
        breadcrumbs={
          <Breadcrumbs items={[{ label: 'Genomes', to: '/genomes' }, { label: 'Build asset' }]} />
        }
        lede={
          <>
            Pick a recipe and a genome, fill in the inputs that recipe asks for, and refgenie
            runs the build and registers what comes out as a new asset. Progress shows up in the
            job console at the bottom of the window.
          </>
        }
      />

      <form
        className="flex flex-col gap-6 max-w-screen-md"
        onSubmit={(event) => void handleSubmit(event)}
      >
        <GenomeSelect
          value={values.genome}
          error={show('genome')}
          onChange={(alias, digest) =>
            setValues((current) => ({ ...current, genome: alias, genomeDigest: digest }))
          }
          onInitGenome={canInitGenome ? () => setInitOpen(true) : undefined}
        />

        <RecipeSelect
          recipeName={values.recipeName}
          recipeVersion={values.recipeVersion}
          onChange={handleRecipeChange}
          error={show('recipe')}
        />

        {/* Everything past the recipe is downstream of it: the inputs it asks
            for, the asset it registers, what happens to that asset next. */}
        <div className="rg-build-branch">
          <RecipeInputFields
            recipe={recipe}
            genomeDigest={values.genomeDigest}
            genomeName={values.genome}
            values={values}
            errors={attempted ? errors : serverErrors}
            onChange={setInput}
            buildLinkFor={(assetClass) =>
              `/build?genome=${encodeURIComponent(values.genome)}` +
              `&recipe=${encodeURIComponent(assetClass)}`
            }
          />

          <section className="rg-build-group" aria-label="Output">
            <h2 className="rg-build-group__caption">Output</h2>

            <FormField
              htmlFor="build-asset-group"
              label="Asset group"
              required
              error={show('asset_group')}
              hint="Defaults to the recipe name, matching `refgenie build genome/asset_group`. Clear it to go back to that default."
            >
              <input
                id="build-asset-group"
                className="rg-field__input"
                type="text"
                value={values.assetGroupName}
                onChange={(event) => {
                  groupEdited.current = event.target.value !== '';
                  setField('assetGroupName', event.target.value);
                }}
              />
            </FormField>

            <FormField
              htmlFor="build-asset-name"
              label="Asset name"
              required={nameRequest !== null}
              error={show('asset')}
              hint={
                nameRequest ?? "Leave empty to let the recipe's default_asset template resolve it."
              }
            >
              <input
                id="build-asset-name"
                className="rg-field__input"
                type="text"
                placeholder={
                  nameRequest
                    ? "Enter a name (couldn't read the tool's version)"
                    : "auto — from the recipe's default_asset"
                }
                value={values.assetName}
                onChange={(event) => setField('assetName', event.target.value)}
              />
            </FormField>

            <FormField
              htmlFor="build-asset-description"
              label="Description"
              hint="Recorded only when this build registers a new asset. If the registry path below already exists, refgenie skips the build entirely and keeps the description that asset already has."
            >
              <textarea
                id="build-asset-description"
                className="rg-field__input"
                rows={2}
                value={values.assetDescription}
                onChange={(event) => setField('assetDescription', event.target.value)}
              />
            </FormField>
          </section>

          <section className="rg-build-group" aria-label="Options">
            <h2 className="rg-build-group__caption">Options</h2>

            {stageFolder && (
              // The toggle sits OUTSIDE the <label>: a button inside it would
              // tick the checkbox on every click.
              <div className="rg-check-row">
                <label className="rg-check" htmlFor="build-stage">
                  <input
                    id="build-stage"
                    type="checkbox"
                    checked={values.stage}
                    aria-describedby={stageHelp.hintProps.id}
                    onChange={(event) => setField('stage', event.target.checked)}
                  />
                  <span>Stage the asset after building</span>
                </label>
                <button {...stageHelp.toggleProps}>?</button>
                <p {...stageHelp.hintProps}>
                  Staging is what makes a built asset reachable from anywhere but this
                  machine: it writes a tar of the asset into the stage folder and records its
                  sha256. That archive is what this server offers for download, what{' '}
                  <code className="rg-code rg-code--inline">refgenie push</code> uploads, and
                  what another machine&rsquo;s{' '}
                  <code className="rg-code rg-code--inline">refgenie pull</code> fetches. It
                  is a separate phase after the build — silent, and slow on large assets.
                </p>
              </div>
            )}

            <label className="rg-check" htmlFor="build-pull-parents">
              <input
                id="build-pull-parents"
                type="checkbox"
                checked={values.pullParents}
                onChange={(event) => setField('pullParents', event.target.checked)}
              />
              <span>Pull missing input assets from a subscribed server instead of failing</span>
            </label>
          </section>
        </div>

        {formLevelErrors.length > 0 && (
          <ul className="rg-list" role="alert">
            {formLevelErrors.map(([field, message]) => (
              <li className="rg-field__error" key={field}>
                {message}
              </li>
            ))}
          </ul>
        )}

        <div className="rg-form-actions">
          <button
            type="submit"
            className="rg-btn rg-btn--primary"
            disabled={!valid || submitting}
          >
            {submitting ? 'Submitting…' : 'Build asset'}
          </button>
          <button type="button" className="rg-btn" onClick={() => navigate(-1)}>
            Cancel
          </button>
        </div>

        {/* Under the button, not above the options: it is what pressing the
            button will do, so it belongs where the eye already is. */}
        <p className="rg-registry-preview">
          <span className="rg-muted">Will create </span>
          <code className="rg-code rg-code--inline">
            {resolvedPath ?? registryPathPreview(values)}
          </code>
        </p>
      </form>

      <InitGenomeModal isOpen={initOpen} onClose={() => setInitOpen(false)} />
    </div>
  );
}
