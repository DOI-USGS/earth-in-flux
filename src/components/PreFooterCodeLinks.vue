<template>
  <div id="code-repository-link-container">
    <a
      :href=gitHubRepositoryLink
      target="_blank"
      aria-label="github link"
    >See the code behind this visualization
      <font-awesome-icon :icon="{ prefix: 'fab', iconName: 'github' }" />
    </a>
  </div>
</template>

<script setup>
  // define props
  defineProps({
    gitHubRepositoryLink:{
      type: String,
      default: 'https://github.com/DOI-USGS'
    },
  })
</script>

<style scoped lang="scss">

  #code-repository-link-container {
    display: flex;
    justify-content: flex-end;
    width: 100%;
    background-color: var(--color-code-links-background);
    margin: 0 auto;
    padding: 0.4rem 1.5rem;
    a {
      color: var(--color-text);
      margin-left: 10px;
      text-decoration: none;
      font-weight: 300;
    }
  }

</style>
