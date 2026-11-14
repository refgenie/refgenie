import { describe, expect, it, vi } from 'vitest';
import { render, screen } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { Tabs } from './Tabs';

const tabs = [
  { id: 'a', label: 'A' },
  { id: 'b', label: 'B', disabled: true },
  { id: 'c', label: 'C', count: 3 },
];

describe('Tabs', () => {
  it('links tabs to panels and marks the selected one', () => {
    render(
      <Tabs idPrefix="t" label="Test" tabs={tabs} selected="c" onSelect={() => {}} />,
    );
    const c = screen.getByRole('tab', { name: 'C 3' });
    expect(c).toHaveAttribute('aria-selected', 'true');
    expect(c).toHaveAttribute('aria-controls', 't-panel-c');
    expect(c).toHaveAttribute('tabindex', '0');
    expect(screen.getByRole('tab', { name: 'A' })).toHaveAttribute('tabindex', '-1');
  });

  it('skips disabled tabs with the arrow keys and wraps', async () => {
    const onSelect = vi.fn();
    render(
      <Tabs idPrefix="t" label="Test" tabs={tabs} selected="a" onSelect={onSelect} />,
    );
    screen.getByRole('tab', { name: 'A' }).focus();
    await userEvent.keyboard('{ArrowRight}');
    expect(onSelect).toHaveBeenLastCalledWith('c');
    await userEvent.keyboard('{ArrowRight}');
    expect(onSelect).toHaveBeenLastCalledWith('a');
    await userEvent.keyboard('{End}');
    expect(onSelect).toHaveBeenLastCalledWith('c');
  });
});
