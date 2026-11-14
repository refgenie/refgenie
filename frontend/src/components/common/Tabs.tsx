/** Accessible tabs: a tablist with arrow-key navigation, plus matching panels. */

import { useRef } from 'react';
import type { KeyboardEvent, ReactNode } from 'react';
import { cn } from '../../utils/cn';

export interface TabItem<T extends string> {
  id: T;
  label: string;
  /** Shown after the label. `null` reserves the space while it loads. */
  count?: number | null;
  disabled?: boolean;
  /** Tooltip, e.g. why a tab is disabled. */
  title?: string;
}

export interface TabsProps<T extends string> {
  /** Prefix for tab and panel element ids. */
  idPrefix: string;
  label: string;
  tabs: TabItem<T>[];
  /** `null` while the default tab is still being worked out. */
  selected: T | null;
  onSelect: (id: T) => void;
  className?: string;
}

function tabId(prefix: string, id: string): string {
  return `${prefix}-tab-${id}`;
}

function panelId(prefix: string, id: string): string {
  return `${prefix}-panel-${id}`;
}

export function Tabs<T extends string>({
  idPrefix,
  label,
  tabs,
  selected,
  onSelect,
  className,
}: TabsProps<T>) {
  const refs = useRef(new Map<string, HTMLButtonElement>());
  const enabled = tabs.filter((tab) => !tab.disabled);
  // One tab is in the page tab order: the selected one, else the first enabled.
  const focusable = selected ?? enabled[0]?.id;

  const onKeyDown = (event: KeyboardEvent<HTMLDivElement>) => {
    if (enabled.length === 0) return;
    const current = enabled.findIndex(
      (tab) => refs.current.get(tab.id) === document.activeElement,
    );
    let next: number | null = null;
    if (event.key === 'ArrowRight') next = (current + 1) % enabled.length;
    else if (event.key === 'ArrowLeft')
      next = (current - 1 + enabled.length) % enabled.length;
    else if (event.key === 'Home') next = 0;
    else if (event.key === 'End') next = enabled.length - 1;
    if (next === null) return;
    event.preventDefault();
    const target = enabled[next];
    refs.current.get(target.id)?.focus();
    onSelect(target.id);
  };

  return (
    <div role="tablist" aria-label={label} className={cn('rg-tabs', className)}
      onKeyDown={onKeyDown}>
      {tabs.map((tab) => {
        const isSelected = tab.id === selected;
        return (
          <button
            key={tab.id}
            ref={(node) => {
              if (node) refs.current.set(tab.id, node);
              else refs.current.delete(tab.id);
            }}
            type="button"
            role="tab"
            id={tabId(idPrefix, tab.id)}
            aria-controls={panelId(idPrefix, tab.id)}
            aria-selected={isSelected}
            tabIndex={tab.id === focusable ? 0 : -1}
            disabled={tab.disabled}
            title={tab.title}
            className={cn('rg-tabs__tab', isSelected && 'rg-tabs__tab--selected')}
            onClick={() => onSelect(tab.id)}
          >
            {tab.label}
            {tab.count !== undefined && (
              <span className="rg-tabs__count" aria-hidden={tab.count === null}>
                {tab.count === null ? '' : tab.count.toLocaleString()}
              </span>
            )}
          </button>
        );
      })}
    </div>
  );
}

export interface TabPanelProps {
  idPrefix: string;
  id: string;
  children: ReactNode;
  className?: string;
}

/** The panel for one tab. Render only the selected one. */
export function TabPanel({ idPrefix, id, children, className }: TabPanelProps) {
  return (
    <div
      role="tabpanel"
      id={panelId(idPrefix, id)}
      aria-labelledby={tabId(idPrefix, id)}
      tabIndex={0}
      className={cn('rg-tabs__panel', className)}
    >
      {children}
    </div>
  );
}
