/**
 * The recipe-driven half of the build form, generated from the three JSON
 * columns on `RecipePublic`.
 *
 * Field names are dotted (`params.x`, `files.x`, `assets.x`) so a server-side
 * `field` on the error envelope lands under the right input with no mapping
 * table in between.
 *
 * `input_files` are **server-side filesystem paths**, not uploads. Local mode
 * means the server and the browser are the same machine, so a path is the
 * correct input; an upload would push a 3 GB FASTA through the HTTP layer to
 * land it on the disk it is already on.
 */

import { useMemo } from 'react';
import { Link } from 'react-router-dom';
import { useAssets } from '../../hooks/queries/useAssets';
import { FormField } from '../common/FormField';
import { inputFieldKey } from './validate';
import type { FieldErrors } from './validate';
import type { RecipePublic } from '../../types/api';

export type InputSection = 'params' | 'files' | 'assets';

export interface RecipeInputFieldsProps {
  recipe: RecipePublic | null;
  genomeDigest: string | null;
  genomeName: string;
  values: {
    params: Record<string, string>;
    files: Record<string, string>;
    assets: Record<string, string>;
  };
  errors: FieldErrors;
  onChange: (section: InputSection, name: string, value: string) => void;
  /** Opens the build form for the missing input's own recipe. */
  buildLinkFor?: (assetClass: string) => string;
}

interface EntrySpec {
  description?: string;
  default?: unknown;
  asset_class?: string;
  [key: string]: unknown;
}

function specEntries(
  source: Record<string, Record<string, unknown>> | null | undefined,
): Array<[string, EntrySpec]> {
  return Object.entries(source ?? {}) as Array<[string, EntrySpec]>;
}

function describe(spec: EntrySpec): string | undefined {
  return typeof spec.description === 'string' ? spec.description : undefined;
}

function defaultText(spec: EntrySpec): string {
  const value = spec.default;
  if (value === undefined || value === null) return '';
  return String(value);
}

export function RecipeInputFields({
  recipe,
  genomeDigest,
  genomeName,
  values,
  errors,
  onChange,
  buildLinkFor,
}: RecipeInputFieldsProps) {
  const inputAssets = specEntries(recipe?.input_assets);

  // The endpoint has no asset-class filter, so the class match is client-side.
  const assets = useAssets(
    { genome_digest: genomeDigest ?? undefined, limit: 200 },
    { enabled: !!genomeDigest && inputAssets.length > 0 },
  );

  const candidates = useMemo(() => {
    const byClass = new Map<string, Array<{ label: string; value: string }>>();
    for (const asset of assets.data?.items ?? []) {
      const className = asset.asset_class_name;
      if (!className || !asset.asset_group_name) continue;
      const list = byClass.get(className) ?? [];
      list.push({
        label: `${asset.asset_group_name}:${asset.name}`,
        // The unambiguous form `BuildParams.assets` documents.
        value: `${genomeName || genomeDigest}/${asset.asset_group_name}:${asset.name}`,
      });
      byClass.set(className, list);
    }
    return byClass;
  }, [assets.data, genomeName, genomeDigest]);

  if (!recipe) return null;

  const params = specEntries(recipe.input_params);
  const files = specEntries(recipe.input_files);

  return (
    <>
      {params.length > 0 && (
        <section className="flex flex-col gap-4">
          <h2 className="text-lg font-semibold">Parameters</h2>
          {params.map(([name, spec]) => {
            const fallback = defaultText(spec);
            const numeric = fallback !== '' && Number.isFinite(Number(fallback));
            return (
              <FormField
                key={name}
                htmlFor={`build-param-${name}`}
                label={name}
                required={fallback === ''}
                hint={describe(spec)}
                error={errors[inputFieldKey('params', name)]}
              >
                <input
                  id={`build-param-${name}`}
                  className="rg-field__input"
                  type={numeric ? 'number' : 'text'}
                  value={values.params[name] ?? fallback}
                  placeholder={fallback || undefined}
                  onChange={(event) => onChange('params', name, event.target.value)}
                />
              </FormField>
            );
          })}
        </section>
      )}

      {files.length > 0 && (
        <section className="flex flex-col gap-4">
          <h2 className="text-lg font-semibold">Files</h2>
          <p className="rg-muted text-sm">
            Paths on the machine running refgenie, not uploads. Existence is checked by
            preflight.
          </p>
          {files.map(([name, spec]) => (
            <FormField
              key={name}
              htmlFor={`build-file-${name}`}
              label={name}
              required={defaultText(spec) === ''}
              hint={describe(spec)}
              error={errors[inputFieldKey('files', name)]}
            >
              <input
                id={`build-file-${name}`}
                className="rg-field__input rg-field__input--mono"
                type="text"
                spellCheck={false}
                placeholder="/absolute/path/to/file"
                value={values.files[name] ?? defaultText(spec)}
                onChange={(event) => onChange('files', name, event.target.value)}
              />
            </FormField>
          ))}
        </section>
      )}

      {inputAssets.length > 0 && (
        <section className="flex flex-col gap-4">
          <h2 className="text-lg font-semibold">Input assets</h2>
          {inputAssets.map(([name, spec]) => {
            const assetClass = typeof spec.asset_class === 'string' ? spec.asset_class : '';
            const options = candidates.get(assetClass) ?? [];
            const preselect =
              options.find((option) => option.label === defaultText(spec))?.value ?? '';

            if (!assets.isPending && options.length === 0) {
              // Chained builds are the normal case; a bare empty dropdown here
              // is a dead end, so name the missing class and offer the build.
              return (
                <div className="rg-field" key={name}>
                  <p className="rg-field__label">{name}</p>
                  <p className="rg-muted text-sm">
                    No <code className="rg-code rg-code--inline">{assetClass}</code> asset for
                    this genome.
                  </p>
                  {buildLinkFor && assetClass && (
                    <Link className="rg-link text-sm" to={buildLinkFor(assetClass)}>
                      Build it first
                    </Link>
                  )}
                  {errors[inputFieldKey('assets', name)] && (
                    <p className="rg-field__error" role="alert">
                      {errors[inputFieldKey('assets', name)]}
                    </p>
                  )}
                </div>
              );
            }

            return (
              <FormField
                key={name}
                htmlFor={`build-asset-${name}`}
                label={name}
                required={defaultText(spec) === ''}
                hint={describe(spec) ?? `Any ${assetClass} asset of this genome.`}
                error={errors[inputFieldKey('assets', name)]}
              >
                <select
                  id={`build-asset-${name}`}
                  className="rg-field__input"
                  value={values.assets[name] ?? preselect}
                  onChange={(event) => onChange('assets', name, event.target.value)}
                >
                  <option value="">
                    {assets.isPending ? 'Loading assets…' : 'Choose an asset…'}
                  </option>
                  {options.map((option) => (
                    <option key={option.value} value={option.value}>
                      {option.label}
                    </option>
                  ))}
                </select>
              </FormField>
            );
          })}
        </section>
      )}
    </>
  );
}
