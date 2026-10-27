import { Link } from 'react-router-dom';
import { useAssetClassIndex } from '../../hooks/queries/useAssetClasses';
import { specEntries } from '../../utils/recipes';
import { DescriptionList } from '../common/DescriptionList';
import { JsonBlock } from '../common/JsonBlock';
import type { AssetClassIndex } from '../../hooks/queries/useAssetClasses';
import type { EntrySpec } from '../../utils/recipes';
import type { InputEntities, RecipePublic } from '../../types/api';

export interface RecipeDetailPanelProps {
  recipe: RecipePublic;
}

function InputEntityList({ entities }: { entities: InputEntities | null }) {
  if (!entities || Object.keys(entities).length === 0) return null;
  return (
    <ul className="flex flex-col gap-1">
      {Object.entries(entities).map(([name, spec]) => (
        <li className="text-sm" key={name}>
          <code className="rg-code rg-code--inline">{name}</code>
          {typeof spec?.description === 'string' && (
            <span className="rg-muted"> — {spec.description}</span>
          )}
          {spec?.default !== undefined && spec.default !== null && (
            <span className="rg-muted"> (default: {String(spec.default)})</span>
          )}
        </li>
      ))}
    </ul>
  );
}

function InputAssetList({
  entities,
  classIndex,
}: {
  entities: InputEntities | null;
  classIndex: AssetClassIndex | undefined;
}) {
  const entries = specEntries(entities);
  if (entries.length === 0) return null;
  return (
    <ul className="flex flex-col gap-1">
      {entries.map(([slot, spec]: [string, EntrySpec]) => {
        // `input_assets` records the class by NAME only, so this resolves through
        // `byName`, which takes the first hit when versions share a name. Nothing
        // better is reachable from the payload.
        const assetClass = spec.asset_class ? classIndex?.byName.get(spec.asset_class) : undefined;
        return (
          <li className="text-sm" key={slot}>
            <code className="rg-code rg-code--inline">{slot}</code>
            <span className="rg-muted"> → </span>
            {assetClass?.id !== null && assetClass?.id !== undefined ? (
              <Link className="rg-link" to={`/asset-classes/${assetClass.id}`}>
                {assetClass.name}
              </Link>
            ) : (
              <span className="rg-muted">{spec.asset_class ?? 'unknown'}</span>
            )}
            {typeof spec.description === 'string' && (
              <span className="rg-muted"> — {spec.description}</span>
            )}
          </li>
        );
      })}
    </ul>
  );
}

export function RecipeDetailPanel({ recipe }: RecipeDetailPanelProps) {
  const classIndex = useAssetClassIndex();

  return (
    <div className="flex flex-col gap-6">
      <DescriptionList
        items={[
          { term: 'Name', value: recipe.name },
          { term: 'Version', value: recipe.version },
          { term: 'Description', value: recipe.description ?? '' },
          {
            term: 'Output asset class',
            value: (() => {
              const outputClass = classIndex.data?.byId.get(recipe.output_asset_class_id);
              return (
                <Link className="rg-link" to={`/asset-classes/${recipe.output_asset_class_id}`}>
                  {/* `RecipePublic` carries only the FK, so the name is joined
                      client-side. The raw id is the fallback, not the label. */}
                  {outputClass
                    ? `${outputClass.name} v${outputClass.version}`
                    : `#${recipe.output_asset_class_id}`}
                </Link>
              );
            })(),
          },
          { term: 'Default asset', value: recipe.default_asset },
          { term: 'Docker image', value: recipe.docker_image ?? '' },
          {
            term: 'Inherent',
            value: recipe.inherent?.length ? recipe.inherent.join(', ') : '',
          },
          {
            term: 'Custom seek keys',
            value:
              recipe.custom_seek_keys && Object.keys(recipe.custom_seek_keys).length > 0 ? (
                <JsonBlock value={recipe.custom_seek_keys} label="Show" />
              ) : (
                ''
              ),
          },
        ]}
      />

      <section>
        <h2 className="text-lg font-semibold mb-2">Command templates</h2>
        {recipe.command_templates.length === 0 ? (
          <p className="rg-muted text-sm">No command templates.</p>
        ) : (
          <div className="flex flex-col gap-2">
            {recipe.command_templates.map((template, index) => (
              <pre className="rg-code" key={index}>
                {template}
              </pre>
            ))}
          </div>
        )}
      </section>

      <section className="grid grid-cols-1 md:grid-cols-3 gap-6">
        <div>
          <h2 className="text-sm font-semibold mb-2">Input params</h2>
          <InputEntityList entities={recipe.input_params} />
        </div>
        <div>
          <h2 className="text-sm font-semibold mb-2">Input files</h2>
          <InputEntityList entities={recipe.input_files} />
        </div>
        <div>
          <h2 className="text-sm font-semibold mb-2">Input assets</h2>
          <InputAssetList entities={recipe.input_assets} classIndex={classIndex.data} />
        </div>
      </section>
    </div>
  );
}
