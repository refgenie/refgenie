import { describe, expect, it } from 'vitest';
import { render, screen } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { FormField } from './FormField';

function renderField(props: Partial<Parameters<typeof FormField>[0]> = {}) {
  return render(
    <FormField htmlFor="demo" label="Genome" hint="Type an alias." {...props}>
      <input id="demo" className="rg-field__input" type="text" />
      <datalist id="demo-options" />
    </FormField>,
  );
}

describe('FormField', () => {
  it('labels the control and describes it with the hint, collapsed or not', async () => {
    renderField();
    const control = screen.getByRole('textbox', { name: /Genome/ });

    // Collapsing the hint is a visual state: it stays in the tree and stays
    // wired to the control, or a screen-reader user loses the explainer.
    const hint = screen.getByText('Type an alias.');
    expect(hint).toHaveAttribute('id', 'demo-hint');
    expect(control).toHaveAttribute('aria-describedby', 'demo-hint');
    expect(hint).toHaveClass('sr-only');

    await userEvent.click(screen.getByRole('button', { name: 'Help with Genome' }));
    expect(screen.getByText('Type an alias.')).not.toHaveClass('sr-only');
    expect(control).toHaveAttribute('aria-describedby', 'demo-hint');
  });

  it('gives the help toggle a name, a target and an expanded state', async () => {
    renderField();
    const toggle = screen.getByRole('button', { name: 'Help with Genome' });

    expect(toggle).toHaveAttribute('aria-controls', 'demo-hint');
    expect(toggle).toHaveAttribute('aria-expanded', 'false');

    // Keyboard, not just pointer: a hover-only tooltip is unreachable here.
    toggle.focus();
    await userEvent.keyboard('{Enter}');
    expect(toggle).toHaveAttribute('aria-expanded', 'true');
    await userEvent.keyboard(' ');
    expect(toggle).toHaveAttribute('aria-expanded', 'false');
  });

  it('describes only the control, not its sibling elements', () => {
    renderField();
    expect(document.getElementById('demo-options')).not.toHaveAttribute('aria-describedby');
  });

  it('offers no toggle when there is no hint', () => {
    renderField({ hint: undefined });
    expect(screen.queryByRole('button')).not.toBeInTheDocument();
    expect(screen.getByRole('textbox')).not.toHaveAttribute('aria-describedby');
  });

  it('keeps the error visible and announced, never behind the toggle', () => {
    renderField({ error: 'No such genome.' });
    const error = screen.getByRole('alert');
    expect(error).toHaveTextContent('No such genome.');
    expect(error).not.toHaveClass('sr-only');
  });
});
