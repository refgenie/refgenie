import type { ApiClient, RequestInitLite } from '../http';
import type {
  RelationshipsExpandedResponse,
  RelationshipsResponse,
} from '../../types/api';

export const getRelationships = (
  c: ApiClient,
  digest: string,
  init?: RequestInitLite,
) =>
  c.get<RelationshipsResponse>(
    `/relationships/${encodeURIComponent(digest)}`,
    { expand: false },
    init,
  );

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
