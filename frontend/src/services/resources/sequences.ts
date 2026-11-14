import type { ApiClient, RequestInitLite } from '../http';
import type { GenomeSequencesResponse, SeqcolComparison } from '../../types/api';

const enc = encodeURIComponent;

export const getGenomeSequences = (c: ApiClient, digest: string, init?: RequestInitLite) =>
  c.get<GenomeSequencesResponse>(`/genomes/${enc(digest)}/sequences`, undefined, init);

export const compareGenomes = (c: ApiClient, a: string, b: string, init?: RequestInitLite) =>
  c.get<SeqcolComparison>(`/genomes/${enc(a)}/compare/${enc(b)}`, undefined, init);
