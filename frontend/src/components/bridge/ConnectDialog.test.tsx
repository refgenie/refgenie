import { beforeEach, describe, expect, it } from 'vitest';
import { screen } from '@testing-library/react';
import { ConnectDialog } from './ConnectDialog';
import { renderWithProviders, serverConfig } from '../../test/renderWithProviders';
import { useBridgeStore } from '../../stores/bridgeStore';

beforeEach(() => {
  localStorage.clear();
  useBridgeStore.getState().reset();
});

describe('ConnectDialog', () => {
  it('says refgenie dash needs the dash extra, for the server version', () => {
    renderWithProviders(<ConnectDialog open onClose={() => {}} />, {
      config: { ...serverConfig, refgenie_version: '1.0.0a2' },
    });
    expect(screen.getByText('pip install "refgenie[dash]>=1.0.0a2"')).toBeInTheDocument();
    expect(screen.getByRole('button', { name: 'Copy install command' })).toBeInTheDocument();
  });
});
