/**
 * BaseModal Component
 *
 * A reusable modal component with built-in accessibility features:
 * - Escape key to close
 * - Click outside to close
 * - Body scroll lock when open
 * - Focus trap (focuses first focusable element)
 *
 * Required CSS: Import modal.css or add modal classes to components.css
 *
 * @example
 * ```tsx
 * <BaseModal
 *   isOpen={showModal}
 *   onClose={() => setShowModal(false)}
 *   title="Edit Item"
 *   size="md"
 * >
 *   <form onSubmit={handleSubmit}>
 *     <div className="form-group">
 *       <label className="form-label">Name</label>
 *       <input className="form-input" type="text" />
 *     </div>
 *     <BaseModal.Footer>
 *       <button type="button" className="btn btn--secondary" onClick={onClose}>
 *         Cancel
 *       </button>
 *       <button type="submit" className="btn btn--primary">
 *         Save
 *       </button>
 *     </BaseModal.Footer>
 *   </form>
 * </BaseModal>
 * ```
 */

import { useEffect, useRef } from 'react';
import type { ReactNode } from 'react';

interface BaseModalProps {
  /** Whether the modal is visible */
  isOpen: boolean;
  /** Callback when modal should close */
  onClose: () => void;
  /** Modal title displayed in header */
  title: string;
  /** Modal width: sm (400px), md (560px), lg (720px) */
  size?: 'sm' | 'md' | 'lg';
  /** Modal content */
  children: ReactNode;
}

export function BaseModal({
  isOpen,
  onClose,
  title,
  size = 'md',
  children
}: BaseModalProps) {
  const modalRef = useRef<HTMLDivElement>(null);

  // Handle escape key
  useEffect(() => {
    if (!isOpen) return;

    const handleEscape = (e: KeyboardEvent) => {
      if (e.key === 'Escape') {
        onClose();
      }
    };

    document.addEventListener('keydown', handleEscape);
    return () => document.removeEventListener('keydown', handleEscape);
  }, [isOpen, onClose]);

  // Prevent body scroll when modal is open
  useEffect(() => {
    if (isOpen) {
      document.body.style.overflow = 'hidden';
    } else {
      document.body.style.overflow = '';
    }
    return () => {
      document.body.style.overflow = '';
    };
  }, [isOpen]);

  // Focus first focusable element when opened
  useEffect(() => {
    if (isOpen && modalRef.current) {
      const focusable = modalRef.current.querySelector<HTMLElement>(
        'input, textarea, select, button'
      );
      focusable?.focus();
    }
  }, [isOpen]);

  if (!isOpen) return null;

  return (
    <div className="modal-overlay" onClick={onClose}>
      <div
        ref={modalRef}
        className={`modal modal--${size}`}
        onClick={(e) => e.stopPropagation()}
      >
        <div className="modal__header">
          <h2 className="modal__title">{title}</h2>
          <button
            type="button"
            className="modal__close"
            onClick={onClose}
            aria-label="Close"
          >
            &times;
          </button>
        </div>
        <div className="modal__body">
          {children}
        </div>
      </div>
    </div>
  );
}

/** Footer subcomponent for consistent button layout */
BaseModal.Footer = function Footer({ children }: { children: ReactNode }) {
  return <div className="modal__footer">{children}</div>;
};

export default BaseModal;
