export interface JsonBlockProps {
  value: unknown;
  label?: string;
}

/** Raw JSON viewer, collapsed by default. No dangerouslySetInnerHTML anywhere. */
export function JsonBlock({ value, label = 'Raw JSON' }: JsonBlockProps) {
  return (
    <details>
      <summary className="cursor-pointer text-sm rg-link">{label}</summary>
      <pre className="rg-code mt-2">{JSON.stringify(value, null, 2)}</pre>
    </details>
  );
}
