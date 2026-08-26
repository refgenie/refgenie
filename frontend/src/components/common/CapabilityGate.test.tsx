import { describe, expect, it } from 'vitest';
import { screen } from '@testing-library/react';
import { CapabilityGate } from './CapabilityGate';
import { localConfig, renderWithProviders, serverConfig } from '../../test/renderWithProviders';

describe('CapabilityGate', () => {
  it('renders children when the capability is on', () => {
    renderWithProviders(
      <CapabilityGate cap="pull">
        <span>pull ui</span>
      </CapabilityGate>,
      { config: localConfig },
    );
    expect(screen.getByText('pull ui')).toBeInTheDocument();
  });

  it('renders nothing when the capability is off', () => {
    renderWithProviders(
      <CapabilityGate cap="pull">
        <span>pull ui</span>
      </CapabilityGate>,
      { config: serverConfig },
    );
    expect(screen.queryByText('pull ui')).not.toBeInTheDocument();
  });

  it('gates on the flag, not the mode: server mode still allows downloads', () => {
    renderWithProviders(
      <CapabilityGate cap="downloads">
        <span>download link</span>
      </CapabilityGate>,
      { config: serverConfig },
    );
    expect(screen.getByText('download link')).toBeInTheDocument();
  });
});
