import './styles/main.css';

import { StrictMode } from 'react';
import { createRoot } from 'react-dom/client';
import { App } from './App';
import { ConfigProvider } from './app/ConfigProvider';
import { ToastProvider } from './components/common/ToastProvider';
import { createAppResourceCache } from './app/resourceCache';
import { ResourceCacheContext } from './services/resourceCacheContext';

const root = document.getElementById('root');
if (!root) throw new Error('missing #root element');

const resourceCache = createAppResourceCache();

createRoot(root).render(
  <StrictMode>
    <ResourceCacheContext.Provider value={resourceCache}>
      <ConfigProvider>
        <ToastProvider>
          <App />
        </ToastProvider>
      </ConfigProvider>
    </ResourceCacheContext.Provider>
  </StrictMode>,
);
