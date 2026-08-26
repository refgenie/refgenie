import { cn } from '../../utils/cn';

export type BadgeVariant = 'file' | 'archive' | 'default' | 'local' | 'server';

export interface BadgeProps {
  variant?: BadgeVariant;
  children: React.ReactNode;
  title?: string;
}

export function Badge({ variant = 'default', children, title }: BadgeProps) {
  return (
    <span className={cn('rg-badge', `rg-badge--${variant}`)} title={title}>
      {children}
    </span>
  );
}
