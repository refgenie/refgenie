import { listQuery } from '../../utils/pagination';
import type { ApiClient, RequestInitLite } from '../http';
import type { ConfigurationPublic } from '../../types/api';
import type { ListParams, PageResponse } from '../../types/pagination';

export const listConfigurations = (
  c: ApiClient,
  p: ListParams = {},
  init?: RequestInitLite,
) => c.get<PageResponse<ConfigurationPublic>>('/configurations', listQuery(p), init);
