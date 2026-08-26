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
import { FormField } from '../components/common/FormField';
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
  const [resolvedPath, setResolvedPath] = useState<string | null>(null);
  const [submitting, setSubmitting] = useState(false);
  const [attempted, setAttempted] = useState(false);
  const [initOpen, setInitOpen] = useState(false);
  const prefilled = useRef(false);

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
    () => validateBuildForm(recipe, values),
    [recipe, values],
  );
  const errors: FieldErrors = { ...serverErrors, ...clientErrors };
  const valid = Object.keys(clientErrors).length === 0;

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
    if (!valid || !canBuild) return;
    const timer = setTimeout(() => {
      preflightBuild(client, toRequest(values))
        .then((response) => {
          const next: FieldErrors = {};
          for (const item of response.errors ?? []) {
            if (item.field) next[item.field] = item.message;
          }
          setServerErrors(next);
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
        });
    }, PREFLIGHT_DEBOUNCE_MS);
    return () => clearTimeout(timer);
  }, [client, values, valid, canBuild]);

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
      assetGroupName: current.assetGroupName || (next?.name ?? ''),
      params: {},
      files: {},
      assets: {},
    }));
    setServerErrors({});
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
    <div className="flex flex-col gap-6 max-w-screen-md">
      <Breadcrumbs items={[{ label: 'Genomes', to: '/genomes' }, { label: 'Build asset' }]} />
      <h1 className="text-3xl font-bold">Build asset</h1>

      <form className="flex flex-col gap-6" onSubmit={(event) => void handleSubmit(event)}>
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

        <section className="flex flex-col gap-4">
          <h2 className="text-lg font-semibold">Output</h2>

          <FormField
            htmlFor="build-asset-group"
            label="Asset group"
            required
            error={show('asset_group')}
            hint="Defaults to the recipe name, matching `refgenie build genome/asset_group`."
          >
            <input
              id="build-asset-group"
              className="rg-field__input"
              type="text"
              value={values.assetGroupName}
              onChange={(event) => setField('assetGroupName', event.target.value)}
            />
          </FormField>

          <FormField
            htmlFor="build-asset-name"
            label="Asset name"
            error={show('asset')}
            hint="Leave empty to let the recipe's default_asset template resolve it."
          >
            <input
              id="build-asset-name"
              className="rg-field__input"
              type="text"
              placeholder="auto — from the recipe's default_asset"
              value={values.assetName}
              onChange={(event) => setField('assetName', event.target.value)}
            />
          </FormField>

          <FormField
            htmlFor="build-asset-description"
            label="Description"
            hint="Only applied when the asset does not already exist."
          >
            <textarea
              id="build-asset-description"
              className="rg-field__input"
              rows={2}
              value={values.assetDescription}
              onChange={(event) => setField('assetDescription', event.target.value)}
            />
          </FormField>

          <p className="rg-registry-preview">
            <span className="rg-muted">Will create </span>
            <code className="rg-code rg-code--inline">
              {resolvedPath ?? registryPathPreview(values)}
            </code>
          </p>
        </section>

        <section className="flex flex-col gap-3">
          <h2 className="text-lg font-semibold">Options</h2>

          {stageFolder && (
            <label className="rg-check" htmlFor="build-stage">
              <input
                id="build-stage"
                type="checkbox"
                checked={values.stage}
                onChange={(event) => setField('stage', event.target.checked)}
              />
              <span>
                Stage the asset after building
                <span className="rg-field__hint">
                  Adds a tar plus a full sha256 of the tarball once the build finishes — a
                  distinct, silent, potentially long phase.
                </span>
              </span>
            </label>
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
            {submitting ? 'Submitting…' : 'Build'}
          </button>
          <button type="button" className="rg-btn" onClick={() => navigate(-1)}>
            Cancel
          </button>
        </div>
      </form>

      <InitGenomeModal isOpen={initOpen} onClose={() => setInitOpen(false)} />
    </div>
  );
}
