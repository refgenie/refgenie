/**
 * The contract the whole link-affordance system rests on.
 *
 * One rule, tested three ways: an icon is never the accessible name of
 * anything. Every glyph is `aria-hidden`, and the fact it carries reaches the
 * accessibility tree as text from the component that owns it.
 */

import { describe, expect, it } from 'vitest';
import { render, screen } from '@testing-library/react';
import { ApiLink } from './ApiLink';
import { CopyButton } from './CopyButton';
import { DownloadLink } from './DownloadLink';
import { ExternalLink } from './ExternalLink';
import { Icon, ICON_NAMES } from './Icon';

describe('Icon', () => {
  it.each(ICON_NAMES)('renders %s as an aria-hidden svg', (name) => {
    const { container } = render(<Icon name={name} />);
    const svg = container.querySelector('svg');
    expect(svg).not.toBeNull();
    expect(svg).toHaveAttribute('aria-hidden', 'true');
    expect(svg).toHaveClass('rg-icon');
    // Colour and size come from CSS, never from the element.
    expect(svg).toHaveAttribute('stroke', 'currentColor');
    expect(svg).not.toHaveAttribute('style');
  });
});

describe('ExternalLink', () => {
  it('announces the new tab in the accessible name, not the glyph', () => {
    render(<ExternalLink href="https://example.org">Documentation</ExternalLink>);
    const link = screen.getByRole('link', { name: 'Documentation (opens in a new tab)' });
    expect(link).toHaveAttribute('target', '_blank');
    expect(link).toHaveAttribute('rel', 'noopener noreferrer');
    expect(link.querySelector('svg')).not.toBeNull();
  });

  it('renders no glyph when the indicator would be decoration', () => {
    const { container } = render(
      <ExternalLink href="https://example.org" icon="none">
        Documentation
      </ExternalLink>,
    );
    expect(container.querySelector('svg')).toBeNull();
  });
});

describe('ApiLink', () => {
  it('opens this instance JSON in a new tab', () => {
    render(<ApiLink href="/openapi.json">OpenAPI schema</ApiLink>);
    const link = screen.getByRole('link', { name: 'OpenAPI schema (opens in a new tab)' });
    expect(link).toHaveAttribute('target', '_blank');
    expect(link).toHaveAttribute('rel', 'noopener noreferrer');
  });
});

describe('DownloadLink', () => {
  it('puts the verb first in the accessible name', () => {
    render(<DownloadLink href="/v4/assets/abc/files/hg38.fa.fai">hg38.fa.fai</DownloadLink>);
    expect(screen.getByRole('link', { name: 'Download hg38.fa.fai' })).toBeInTheDocument();
  });

  it('opens no tab and sets no download attribute', () => {
    // `download` is ignored cross-origin, which is exactly the localhost-bridge
    // case; `target` would flash a tab open and closed on an attachment.
    render(<DownloadLink href="/v4/archives/abc">Archive</DownloadLink>);
    const link = screen.getByRole('link', { name: 'Download Archive' });
    expect(link).not.toHaveAttribute('target');
    expect(link).not.toHaveAttribute('download');
  });
});

describe('CopyButton', () => {
  it('keeps its label after going icon-only', () => {
    render(<CopyButton value="pip install refgenie" label="Copy: pip install refgenie" />);
    const button = screen.getByRole('button', { name: 'Copy: pip install refgenie' });
    expect(button.querySelector('svg')).not.toBeNull();
    expect(button).toHaveTextContent('');
  });
});
