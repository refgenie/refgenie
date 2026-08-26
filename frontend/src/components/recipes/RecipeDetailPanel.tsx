import { Link } from 'react-router-dom';
import { DescriptionList } from '../common/DescriptionList';
import { JsonBlock } from '../common/JsonBlock';
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

export function RecipeDetailPanel({ recipe }: RecipeDetailPanelProps) {
  return (
    <div className="flex flex-col gap-6">
      <DescriptionList
        items={[
          { term: 'Name', value: recipe.name },
          { term: 'Version', value: recipe.version },
          { term: 'Description', value: recipe.description ?? '' },
          {
            term: 'Output asset class',
            value: (
              <Link className="rg-link" to={`/asset-classes/${recipe.output_asset_class_id}`}>
                #{recipe.output_asset_class_id}
              </Link>
            ),
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
          <InputEntityList entities={recipe.input_assets} />
        </div>
      </section>
    </div>
  );
}
