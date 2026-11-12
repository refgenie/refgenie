import { describe, expect, it } from 'vitest';
import { screen } from '@testing-library/react';
import { NoAssetsPanel } from './NoAssetsPanel';
import { renderWithProviders, serverConfig } from '../../test/renderWithProviders';

const D = '5cb5b9d2d1e2b5e5e2f0a3c4d5e6f708';

describe('NoAssetsPanel', () => {
  it('links to the prefilled asset request form', () => {
    renderWithProviders(<NoAssetsPanel digest={D} alias="sacCer3" region="chrI:0-100" />, {
      config: serverConfig,
    });
    expect(screen.getByRole('link', { name: /Request assets/ })).toHaveAttribute(
      'href',
      `https://github.com/refgenie/refgenie-registry/issues/new?template=asset_request.yml&genome=sacCer3&genome_digest=${D}`,
    );
  });

  it('shows both command lines with an alias and a region', () => {
    renderWithProviders(<NoAssetsPanel digest={D} alias="sacCer3" region="chrI:0-100" />, {
      config: serverConfig,
    });
    expect(screen.getByText('refgenie pull --init -g sacCer3')).toBeInTheDocument();
    expect(
      screen.getByText('refgenie getseq -g sacCer3 -l chrI:0-100'),
    ).toBeInTheDocument();
  });

  it('shows no command lines without an alias', () => {
    renderWithProviders(<NoAssetsPanel digest={D} region="chrI:0-100" />, {
      config: serverConfig,
    });
    expect(screen.queryByText(/refgenie pull/)).not.toBeInTheDocument();
    expect(screen.queryByText(/Only need the sequence/)).not.toBeInTheDocument();
  });
});
