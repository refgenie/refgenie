/**
 * A job's log tail. Shows refgenie's own log lines and, for builds, the tool
 * output from the pipeline log, styled differently.
 */

import { useEffect, useRef } from 'react';
import { useJobStore } from '../../stores/jobStore';
import { cn } from '../../utils/cn';
import type { JobStatus } from '../../services/contracts';

export interface JobLogTailProps {
  jobId: string;
  status: JobStatus;
  /** Builds get a taller tail: subprocess output is all they have. */
  lines?: number;
}

export function JobLogTail({ jobId, status, lines = 12 }: JobLogTailProps) {
  // A selector so only this component re-renders as log lines arrive.
  const log = useJobStore((state) => state.logs[jobId]);
  const containerRef = useRef<HTMLDivElement>(null);
  const stickToBottom = useRef(true);

  const tail = log ? log.slice(Math.max(0, log.length - lines)) : [];

  useEffect(() => {
    const node = containerRef.current;
    if (!node || !stickToBottom.current) return;
    node.scrollTop = node.scrollHeight;
  }, [tail.length]);

  // Stop auto-scrolling as soon as the user scrolls up.
  const handleScroll = () => {
    const node = containerRef.current;
    if (!node) return;
    const distance = node.scrollHeight - node.scrollTop - node.clientHeight;
    stickToBottom.current = distance < 8;
  };

  if (status === 'queued') return null;

  if (tail.length === 0) {
    // A build's log starts only once the asset name is resolved.
    return <p className="rg-log-tail__empty rg-muted">Waiting for output…</p>;
  }

  return (
    <div
      className="rg-log-tail"
      ref={containerRef}
      onScroll={handleScroll}
      aria-live="off"
      role="log"
    >
      {tail.map((line, index) => (
        <span
          className={cn('rg-log-tail__line', line.source === 'pipeline' && 'rg-log-tail__line--pipeline')}
          key={`${index}-${line.text}`}
        >
          {line.text}
        </span>
      ))}
    </div>
  );
}
