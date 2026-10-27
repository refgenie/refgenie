import { Badge } from '../common/Badge';
import type { BadgeVariant } from '../common/Badge';

export interface ServingModeBadgesProps {
  modes: string[] | null | undefined;
}

function variantFor(mode: string): BadgeVariant {
  if (mode === 'file') return 'file';
  if (mode === 'archive') return 'archive';
  return 'default';
}

/** The *resolved* serving modes (override → asset class → default). */
export function ServingModeBadges({ modes }: ServingModeBadgesProps) {
  if (!modes || modes.length === 0) return <span className="rg-muted">NA</span>;
  return (
    <span className="flex flex-wrap gap-1">
      {modes.map((mode) => (
        <Badge key={mode} variant={variantFor(mode)}>
          {mode}
        </Badge>
      ))}
    </span>
  );
}
