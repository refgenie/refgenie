import { describe, expect, it } from 'vitest';
import { screen } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { JsonBlock } from './JsonBlock';
import { renderWithProviders } from '../../test/renderWithProviders';

describe('JsonBlock', () => {
  it('opens the JSON in a modal and closes it again', async () => {
    renderWithProviders(<JsonBlock value={{ a: 1 }} label="Show" title="Custom seek keys" />);
    expect(screen.queryByRole('dialog')).toBeNull();

    await userEvent.click(screen.getByRole('button', { name: 'Show' }));
    const dialog = screen.getByRole('dialog', { name: 'Custom seek keys' });
    expect(dialog).toHaveTextContent('"a": 1');

    await userEvent.click(screen.getByRole('button', { name: 'Close' }));
    expect(screen.queryByRole('dialog')).toBeNull();
  });

  it('titles the modal with the label when no title is given', async () => {
    renderWithProviders(<JsonBlock value={[1, 2]} />);
    await userEvent.click(screen.getByRole('button', { name: 'Raw JSON' }));
    expect(screen.getByRole('dialog', { name: 'Raw JSON' })).toBeInTheDocument();
  });
});
