// The types the sequence views consume, kept in one module so a future shared
// package (refget's Explorer and refgenie) can export exactly these.
export type {
  AliasEntry,
  SeqcolComparison,
  SequenceResponse as SequenceRow,
} from '../../types/api';
