/**
 * The progress bar.
 *
 * `percent: null` is the COMMON case, not the exception: builds have no
 * percentage at all, and a pull only has one when the JobManager instrumented
 * the download. So the indeterminate variant is a first-class rendering with
 * the phase name as its label, not a fallback.
 *
 * A native `<progress>` element carries the fill, which is what lets this be
 * both accessible and free of inline styles — a width computed into a `style`
 * attribute is exactly what the style rules forbid.
 */

import { cn } from '../../utils/cn';
import { formatBytes } from '../../utils/format';
import { phaseInfo, phaseStep } from './phases';
import type { JobKind, JobProgressInfo } from '../../services/contracts';

export interface JobProgressProps {
  kind: JobKind;
  progress: JobProgressInfo | null;
  /** Inline variant for the pull button's mini-bar. */
  compact?: boolean;
}

export function JobProgress({ kind, progress, compact = false }: JobProgressProps) {
  const phase = progress?.phase ?? null;
  const info = phaseInfo(kind, phase);
  const step = phaseStep(kind, phase);
  const percent = progress?.percent ?? null;
  const determinate = typeof percent === 'number' && Number.isFinite(percent);
  const rounded = determinate ? Math.max(0, Math.min(100, Math.round(percent))) : null;

  const bytesDone = progress?.bytes_done ?? null;
  const bytesTotal = progress?.bytes_total ?? null;
  const bytesLabel =
    bytesDone === null
      ? null
      : bytesTotal === null
        ? formatBytes(bytesDone)
        : `${formatBytes(bytesDone)} / ${formatBytes(bytesTotal)}`;

  const label = determinate ? `${info.label} — ${rounded}%` : info.label;

  return (
    <div className={cn('rg-progress', compact && 'rg-progress--compact')}>
      {determinate ? (
        <progress
          className="rg-progress__bar"
          value={rounded ?? 0}
          max={100}
          aria-label={label}
          aria-valuenow={rounded ?? 0}
          aria-valuemin={0}
          aria-valuemax={100}
        />
      ) : (
        <div
          className="rg-progress__bar rg-progress__bar--indeterminate"
          role="progressbar"
          aria-label={label}
        >
          <span className="rg-progress__pulse" />
        </div>
      )}

      {!compact && (
        <div className="rg-progress__meta">
          <span className="rg-progress__phase">
            {info.label}
            {step && <span className="rg-muted"> — step {step.index} of {step.total}</span>}
          </span>
          {bytesLabel && <span className="rg-progress__bytes">{bytesLabel}</span>}
        </div>
      )}

      {!compact && info.slowNote && (
        <p className="rg-progress__note rg-muted">{info.slowNote}</p>
      )}
      {!compact && progress?.message && (
        <p className="rg-progress__message rg-muted">{progress.message}</p>
      )}
    </div>
  );
}
