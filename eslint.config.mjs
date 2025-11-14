import pluginVue from 'eslint-plugin-vue'
import { defineConfig, globalIgnores } from "eslint/config";

export default defineConfig( [
  // add rulesets here
  ...pluginVue.configs['flat/recommended'],
  globalIgnores(
    [
      ".gitignore",
      ".vite/"
    ]
  ),
  {
    rules: {
      // disabling v-html warning b/c we provide all injected html
      'vue/no-v-html': 'off'
    }
  }
]);