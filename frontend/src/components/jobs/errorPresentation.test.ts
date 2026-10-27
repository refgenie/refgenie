import { describe, expect, it } from 'vitest';
import { JOB_ERROR_PRESENTATION, presentJobError } from './errorPresentation';

describe('JOB_ERROR_PRESENTATION', () => {
  it('gives every mapped code a headline and a tone', () => {
    for (const [code, presentation] of Object.entries(JOB_ERROR_PRESENTATION)) {
      expect(presentation.headline, code).toBeTruthy();
      expect(presentation.tone, code).toBeTruthy();
    }
  });

  it('offers overwrite for an already-downloaded asset', () => {
    const presentation = presentJobError({ code: 'asset_exists', message: 'exists' });
    expect(presentation.tone).toBe('warning');
    expect(presentation.action?.id).toBe('pull_force');
  });

  it('leads a build failure with the log, since there is no exception behind it', () => {
    const presentation = presentJobError({ code: 'build_failed', message: 'failed' });
    expect(presentation.preformatted).toBe(true);
    expect(presentation.action?.id).toBe('view_log');
  });

  it('explains a digest mismatch instead of offering a retry that would fail again', () => {
    const presentation = presentJobError({ code: 'remote_digest_mismatch', message: 'x' });
    expect(presentation.action).toBeUndefined();
    expect(presentation.note).toBeTruthy();
  });

  it('falls through to the message plus a disclosure for an unknown code', () => {
    const presentation = presentJobError({ code: 'brand_new_code', message: 'Boom' });
    expect(presentation.headline).toBe('Boom');
    expect(presentation.disclosure).toBe(true);
    expect(presentation.action?.id).toBe('retry');
  });
});
