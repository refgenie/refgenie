import type { ApiClient, RequestInitLite } from '../http';
import type { RelationshipsExpandedResponse } from '../../types/api';

export const getRelationshipsExpanded = (
  c: ApiClient,
  digest: string,
  init?: RequestInitLite,
) =>
  c.get<RelationshipsExpandedResponse>(
    `/relationships/${encodeURIComponent(digest)}`,
    { expand: true },
    init,
  );
