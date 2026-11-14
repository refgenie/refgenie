import { useState } from 'react';
import type { ReactNode } from 'react';
import { Copyable } from './Copyable';
import { TabPanel, Tabs } from './Tabs';

export interface CodeSnippet {
  id: string;
  /** The tab label: a language or tool, e.g. "Python". */
  label: string;
  /** How to get the tool, e.g. `pip install refget`. */
  install?: string;
  code: string;
  /** One short line under the code. */
  note?: ReactNode;
}

export interface CodeSnippetsProps {
  /** Prefix for tab ids; unique on the page. */
  idPrefix: string;
  /** Accessible name of the language tabs. */
  label: string;
  snippets: CodeSnippet[];
}

/** One snippet at a time, picked by language, each part with its own copy button. */
export function CodeSnippets({ idPrefix, label, snippets }: CodeSnippetsProps) {
  const [chosen, setChosen] = useState(snippets[0]?.id ?? '');
  const snippet = snippets.find((s) => s.id === chosen) ?? snippets[0];
  if (!snippet) return null;
  return (
    <div className="rg-snippets">
      <Tabs
        idPrefix={idPrefix}
        label={label}
        tabs={snippets.map(({ id, label: tab }) => ({ id, label: tab }))}
        selected={snippet.id}
        onSelect={setChosen}
        className="rg-tabs--sm"
      />
      <TabPanel idPrefix={idPrefix} id={snippet.id} className="rg-snippets__panel">
        {snippet.install && (
          <Copyable
            variant="block"
            value={snippet.install}
            label={`Copy the ${snippet.label} install command`}
          />
        )}
        <Copyable
          variant="block"
          value={snippet.code}
          label={`Copy the ${snippet.label} code`}
        />
        {snippet.note && <p className="rg-snippets__note">{snippet.note}</p>}
      </TabPanel>
    </div>
  );
}
