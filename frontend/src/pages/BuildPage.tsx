/** The build form. A full page, not a modal, because the recipe decides how many fields it has. */

import { useCallback, useEffect, useMemo, useRef, useState } from 'react';
import { useNavigate, useSearchParams } from 'react-router-dom';
import { useLocalApiClient } from '../hooks/useApiClient';
import { useCapability } from '../hooks/useCapability';
import { useConfigurations } from '../hooks/queries/useConfigurations';
import { useToast } from '../hooks/useToast';
import { buildAsset, preflightBuild } from '../services/actions';
import { ApiError } from '../services/http';
import { trackJob } from '../stores/jobStore';
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
import type { RecipeResponse } from '../types/api';

const PREFLIGHT_DEBOUNCE_MS = 400;

/** Form state to request body. The backend forbids extra fields; recipe inputs nest under `params`. */
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
  // Set when preflight says `asset_name_required` (the recipe could not read a
  // tool version to name the asset). Makes Asset name required; holds the
  // server's explanation for the field hint.
  const [nameRequest, setNameRequest] = useState<string | null>(null);
  const [resolvedPath, setResolvedPath] = useState<string | null>(null);
  const [submitting, setSubmitting] = useState(false);
  const [attempted, setAttempted] = useState(false);
  const [initOpen, setInitOpen] = useState(false);
  const prefilled = useRef(false);
  // True once the user types an asset group. Until then (or after clearing it)
  // the group follows the selected recipe.
  const groupEdited = useRef(false);
  const stageHelp = useHelpDisclosure('build-stage-hint', 'staging');

  // Staging fails without a stage folder, so the checkbox only shows when one is set.
  const configurations = useConfigurations({ pageSize: 1 });
  const stageFolder = configurations.data?.results?.[0]?.genome_stage_folder ?? null;

  const recipe: RecipeResponse | null = useMemo(() => {
    if (!values.recipeName) return null;
    const versions = catalog.byName.get(values.recipeName) ?? [];
    return versions.find((r) => r.version === values.recipeVersion) ?? versions[0] ?? null;
  }, [catalog, values.recipeName, values.recipeVersion]);

  // Prefill once, from ?genome=&recipe=&asset_group=, so "build this" links work.
  useEffect(() => {
    if (prefilled.current || catalog.isPending) return;
    prefilled.current = true;
    // `genome` is an alias (preflight does not accept a digest); the digest is
    // only used for the job target key.
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
  // A missing required name blocks submit but not preflight, which is what
  // clears the requirement once a name is typed.
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

  // Debounced server preflight catches missing paths and parent assets.
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
          // Preflight is advisory; the job reports any real failure.
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

  const handleRecipeChange = (next: RecipeResponse | null) => {
    setValues((current) => ({
      ...current,
      recipeName: next?.name ?? null,
      recipeVersion: next?.version ?? null,
      // As in the CLI, the group defaults to the recipe name.
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
      trackJob(
        ref,
        {
          noun: 'build',
          label: `build ${registryPathPreview(values)}`,
          target: {
            genome_digest: values.genomeDigest,
            genome_name: values.genome,
            asset_group_name: values.assetGroupName,
            asset_name: values.assetName || null,
          },
        },
        toast,
      );
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

  // Preflight errors that no input owns (e.g. `stage`, `params`) show above the submit button.
  const claimsAnInput = (field: string) =>
    ['genome', 'recipe', 'asset_group', 'asset'].includes(field) ||
    /^params\.(params|files|assets)\./.test(field);
  const formLevelErrors = Object.entries(serverErrors).filter(
    ([field]) => !claimsAnInput(field),
  );

  return (
    // Cap the form's width, not the MiniHero's parent (see MiniHero).
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
              // The "?" button sits outside the <label> so clicking it does not tick the box.
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
