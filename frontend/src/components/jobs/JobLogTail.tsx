/**
 * The log tail: the primary feedback surface, because there is no progress
 * callback in refgenie for builds and only an optional one for pulls.
 *
 * Two sources reach it and both are legitimate — refgenie's own `logging`
 * lines, demultiplexed per job by the JobManager, and (for builds) the tee'd
 * subprocess output from the pypiper log file. They render in one stream with
 * the pipeline lines visually distinct, so a wall of tool output does not read
 * as refgenie's own phase messages.
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
  // Subscribed via a selector so this component re-renders alone when log
  // lines arrive at tens of events per second.
  const log = useJobStore((state) => state.logs[jobId]);
  const containerRef = useRef<HTMLDivElement>(null);
  const stickToBottom = useRef(true);

  const tail = log ? log.slice(Math.max(0, log.length - lines)) : [];

  useEffect(() => {
    const node = containerRef.current;
    if (!node || !stickToBottom.current) return;
    node.scrollTop = node.scrollHeight;
  }, [tail.length]);

  // Auto-scrolling out from under a reader is the classic bug here, so the
  // stick-to-bottom flag is dropped the moment the user scrolls up.
  const handleScroll = () => {
    const node = containerRef.current;
    if (!node) return;
    const distance = node.scrollHeight - node.scrollTop - node.clientHeight;
    stickToBottom.current = distance < 8;
  };

  if (status === 'queued') return null;

  if (tail.length === 0) {
    // A build's log file does not exist until the asset name is resolved, so
    // the first seconds legitimately have nothing. Never a blank panel.
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
