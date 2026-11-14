import { useEffect, useMemo, useState } from 'react';
import type { KeyboardEvent } from 'react';
import { cn } from '../../utils/cn';

export interface ComboOption {
  value: string;
  label: string;
  /** Other words that find this option, e.g. a genome's other aliases. */
  keywords?: readonly string[];
}

export interface FilterComboboxProps {
  /** The input's id, for the field label. */
  id: string;
  options: readonly ComboOption[];
  /** The chosen option's value. */
  value: string | undefined;
  onChange: (value: string) => void;
  placeholder?: string;
  /** What one option is, for the match count: "genome". */
  noun?: string;
}

interface Match {
  option: ComboOption;
  /** The keyword that matched, when the label did not. */
  via?: string;
}

/** Options whose label or a keyword contains the query, in their given order. */
function filterOptions(options: readonly ComboOption[], query: string): Match[] {
  const needle = query.trim().toLowerCase();
  if (!needle) return options.map((option) => ({ option }));
  const matches: Match[] = [];
  for (const option of options) {
    if (option.label.toLowerCase().includes(needle)) {
      matches.push({ option });
      continue;
    }
    const via = option.keywords?.find((word) => word.toLowerCase().includes(needle));
    if (via) matches.push({ option, via });
  }
  return matches;
}

/**
 * A text box that filters a list of options as you type (the ARIA combobox
 * pattern with list autocomplete). Arrow keys move through the list, Enter
 * picks, Escape closes. Picking calls `onChange`; typing alone does not.
 */
export function FilterCombobox({
  id,
  options,
  value,
  onChange,
  placeholder,
  noun = 'option',
}: FilterComboboxProps) {
  const listId = `${id}-list`;
  const optionId = (index: number) => `${id}-option-${index}`;
  const selectedLabel = options.find((o) => o.value === value)?.label ?? '';

  const [text, setText] = useState(selectedLabel);
  const [query, setQuery] = useState('');
  const [open, setOpen] = useState(false);
  const [active, setActive] = useState(-1);

  // Follow a choice made elsewhere, unless the user is typing.
  const [shownLabel, setShownLabel] = useState(selectedLabel);
  if (selectedLabel !== shownLabel) {
    setShownLabel(selectedLabel);
    if (!open) setText(selectedLabel);
  }

  const matches = useMemo(() => filterOptions(options, query), [options, query]);

  useEffect(() => {
    if (!open || active < 0) return;
    document.getElementById(`${id}-option-${active}`)?.scrollIntoView?.({ block: 'nearest' });
  }, [open, active, id]);

  const openList = () => {
    if (open) return;
    // Opening shows every option, with the current choice active.
    setQuery('');
    setActive(Math.max(0, options.findIndex((o) => o.value === value)));
    setOpen(true);
  };

  const close = () => {
    setOpen(false);
    setActive(-1);
    setText(selectedLabel);
  };

  const pick = (option: ComboOption) => {
    setOpen(false);
    setActive(-1);
    setText(option.label);
    if (option.value !== value) onChange(option.value);
  };

  const onKeyDown = (event: KeyboardEvent<HTMLInputElement>) => {
    const last = matches.length - 1;
    switch (event.key) {
      case 'ArrowDown':
        event.preventDefault();
        if (!open) openList();
        else setActive((i) => (i >= last ? 0 : i + 1));
        break;
      case 'ArrowUp':
        event.preventDefault();
        if (!open) openList();
        else setActive((i) => (i <= 0 ? last : i - 1));
        break;
      case 'Enter':
        if (open && matches[active]) {
          event.preventDefault();
          pick(matches[active].option);
        }
        break;
      case 'Escape':
        if (open) {
          event.preventDefault();
          close();
        }
        break;
      default:
        break;
    }
  };

  return (
    <div className="rg-combo">
      <input
        id={id}
        className="rg-field__input"
        type="text"
        role="combobox"
        aria-autocomplete="list"
        aria-expanded={open}
        aria-controls={listId}
        aria-activedescendant={open && matches[active] ? optionId(active) : undefined}
        autoComplete="off"
        spellCheck={false}
        placeholder={placeholder}
        value={text}
        onChange={(event) => {
          setText(event.target.value);
          setQuery(event.target.value);
          setActive(0);
          setOpen(true);
        }}
        // Selected on focus, so typing replaces the current choice.
        onFocus={(event) => event.target.select()}
        onClick={openList}
        onKeyDown={onKeyDown}
        onBlur={close}
      />
      <ul
        id={listId}
        role="listbox"
        className="rg-combo__list"
        hidden={!open || matches.length === 0}
      >
        {matches.map(({ option, via }, index) => (
          <li
            key={option.value}
            id={optionId(index)}
            role="option"
            aria-selected={index === active}
            className={cn(
              'rg-combo__option',
              index === active && 'rg-combo__option--active',
              option.value === value && 'rg-combo__option--chosen',
            )}
            // Keep focus in the input, so blur does not close before the click.
            onMouseDown={(event) => event.preventDefault()}
            onMouseMove={() => setActive(index)}
            onClick={() => pick(option)}
          >
            {option.label}
            {via && <span className="rg-muted"> · {via}</span>}
          </li>
        ))}
      </ul>
      {open && matches.length === 0 && (
        <p className="rg-combo__empty rg-muted text-sm">No {noun}s match.</p>
      )}
      <span className="sr-only" aria-live="polite">
        {open ? `${matches.length} ${noun}${matches.length === 1 ? '' : 's'}` : ''}
      </span>
    </div>
  );
}
