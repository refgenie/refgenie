import { describe, expect, it } from 'vitest';
import { screen } from '@testing-library/react';
import { MiniHero } from './MiniHero';
import { renderWithProviders } from '../../test/renderWithProviders';

const SERVICE = 'refgenie local dashboard';

describe('MiniHero', () => {
  it('renders the title as the page h1', () => {
    renderWithProviders(<MiniHero title="Recipes" />);
    expect(screen.getByRole('heading', { level: 1, name: 'Recipes' })).toBeInTheDocument();
  });

  it('titles the tab from the heading when no override is given', () => {
    renderWithProviders(<MiniHero title="Recipes" />);
    expect(document.title).toBe(`Recipes · ${SERVICE}`);
  });

  it('honours a documentTitle that differs from the heading', () => {
    renderWithProviders(<MiniHero title="Remote assets" documentTitle="Remote" />);
    expect(screen.getByRole('heading', { level: 1 })).toHaveTextContent('Remote assets');
    expect(document.title).toBe(`Remote · ${SERVICE}`);
  });

  it('falls back to the bare service name while a record is still loading', () => {
    renderWithProviders(<MiniHero title="Genome" documentTitle={false} />);
    expect(document.title).toBe(SERVICE);
  });

  it('renders no lede paragraph and no actions cluster when neither is given', () => {
    const { container } = renderWithProviders(<MiniHero title="Aliases" />);
    expect(container.querySelector('.rg-minihero__lede')).toBeNull();
    expect(container.querySelector('.rg-minihero__actions')).toBeNull();
  });

  it('renders the lede, the actions and the breadcrumb slot when given', () => {
    const { container } = renderWithProviders(
      <MiniHero
        title="Aliases"
        lede={<>An alias is a human-readable name.</>}
        actions={<button type="button">Do a thing</button>}
        breadcrumbs={<nav aria-label="Breadcrumb">crumbs</nav>}
      />,
    );
    expect(container.querySelector('.rg-minihero__lede')).toHaveTextContent(
      'An alias is a human-readable name.',
    );
    expect(screen.getByRole('button', { name: 'Do a thing' })).toBeInTheDocument();
    expect(screen.getByRole('navigation', { name: 'Breadcrumb' })).toBeInTheDocument();
  });
});
