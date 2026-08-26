import './styles/main.css';

import { StrictMode } from 'react';
import { createRoot } from 'react-dom/client';
import { QueryClientProvider } from '@tanstack/react-query';
import { App } from './App';
import { ConfigProvider } from './app/ConfigProvider';
import { ToastProvider } from './components/common/ToastProvider';
import { createQueryClient } from './app/queryClient';

const root = document.getElementById('root');
if (!root) throw new Error('missing #root element');

createRoot(root).render(
  <StrictMode>
    <QueryClientProvider client={createQueryClient()}>
      <ConfigProvider>
        <ToastProvider>
          <App />
        </ToastProvider>
      </ConfigProvider>
    </QueryClientProvider>
  </StrictMode>,
);
