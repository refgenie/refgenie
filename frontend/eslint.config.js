import js from '@eslint/js';
import globals from 'globals';
import reactHooks from 'eslint-plugin-react-hooks';
import reactRefresh from 'eslint-plugin-react-refresh';
import tseslint from 'typescript-eslint';

export default tseslint.config(
  { ignores: ['dist', 'node_modules', 'src/styles/utilities.css', 'src/styles/modal.css'] },
  {
    extends: [js.configs.recommended, ...tseslint.configs.recommended],
    files: ['**/*.{ts,tsx}'],
    languageOptions: {
      ecmaVersion: 2020,
      globals: globals.browser,
    },
    plugins: {
      'react-hooks': reactHooks,
      'react-refresh': reactRefresh,
    },
    rules: {
      ...reactHooks.configs.recommended.rules,
      'react-refresh/only-export-components': ['warn', { allowConstantExport: true }],
      '@typescript-eslint/no-unused-vars': [
        'error',
        { argsIgnorePattern: '^_', varsIgnorePattern: '^_' },
      ],
      'no-restricted-syntax': [
        'error',
        {
          selector: "JSXAttribute[name.name='style']",
          message:
            'No inline styles. Use utility classes or a BEM component in components.css.',
        },
        {
          selector:
            "Program > VariableDeclaration > VariableDeclarator[id.name='API_BASE']",
          message:
            'No module-level API_BASE. The base URL comes from /service-info via ConfigProvider; use useApiClient().',
        },
      ],
    },
  },
  {
    // The route table is JSX that exports a data structure, not a component.
    files: ['src/app/routes.tsx'],
    rules: { 'react-refresh/only-export-components': 'off' },
  },
  {
    // One HTTP layer. A component that calls `fetch` itself bypasses the error
    // envelope, the timeout, and — for mutations — the X-Refgenie-Action header
    // that is the whole cross-origin defence.
    files: ['src/**/*.{ts,tsx}'],
    ignores: ['src/services/**', 'src/test/**'],
    rules: {
      'no-restricted-globals': [
        'error',
        {
          name: 'fetch',
          message:
            'No component calls fetch directly. Use the ApiClient from useApiClient()/useLocalApiClient(), or add a function under src/services/.',
        },
      ],
    },
  },
  {
    files: ['scripts/**/*.mjs'],
    languageOptions: { globals: globals.node },
  },
);
