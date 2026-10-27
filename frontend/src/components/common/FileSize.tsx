import { formatBytes } from '../../utils/format';

export interface FileSizeProps {
  bytes: number | null | undefined;
}

export function FileSize({ bytes }: FileSizeProps) {
  const text = formatBytes(bytes);
  return (
    <span title={bytes === null || bytes === undefined ? undefined : `${bytes} bytes`}>
      {text}
    </span>
  );
}
