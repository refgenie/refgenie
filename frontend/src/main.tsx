import './styles/main.css';

import { StrictMode } from 'react';
import { createRoot } from 'react-dom/client';
import { App } from './App';
import { ConfigProvider } from './app/ConfigProvider';
import { ToastProvider } from './components/common/ToastProvider';
import { appCache } from './app/resourceCache';
import { ResourceCacheContext } from './hooks/useResource';
import { removeLegacyDownloadWorker } from './services/streamDownload';

// An earlier build registered the download worker for the whole site.
removeLegacyDownloadWorker().catch(() => {});

const root = document.getElementById('root');
if (!root) throw new Error('missing #root element');

createRoot(root).render(
  <StrictMode>
    <ResourceCacheContext.Provider value={appCache}>
      <ConfigProvider>
        <ToastProvider>
          <App />
        </ToastProvider>
      </ConfigProvider>
    </ResourceCacheContext.Provider>
  </StrictMode>,
);
