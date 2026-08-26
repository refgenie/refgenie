import { Badge } from '../common/Badge';
import { useUiConfig } from '../../hooks/useUiConfig';

/**
 * The ONLY place the mode string drives the UI (that and the document title).
 * Every affordance is gated on a capability flag instead.
 */
export function ModeBadge() {
  const { mode } = useUiConfig();
  return (
    <Badge variant={mode === 'local' ? 'local' : 'server'} title={`Running in ${mode} mode`}>
      {mode}
    </Badge>
  );
}
