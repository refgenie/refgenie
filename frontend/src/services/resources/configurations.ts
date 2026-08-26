import type { ApiClient, RequestInitLite } from '../http';
import type { ConfigurationPublic } from '../../types/api';
import type { ListParams, Paginated } from '../../types/pagination';

export const listConfigurations = (
  c: ApiClient,
  p: ListParams = {},
  init?: RequestInitLite,
) =>
  c.get<Paginated<ConfigurationPublic>>(
    '/configurations',
    { offset: p.offset, limit: p.limit },
    init,
  );

export const getConfiguration = (c: ApiClient, id: number, init?: RequestInitLite) =>
  c.get<ConfigurationPublic>(`/configurations/${id}`, undefined, init);
