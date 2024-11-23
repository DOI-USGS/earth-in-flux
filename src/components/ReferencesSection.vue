<template>
  <VizSection
    id="references"
    :figures="false"
    :fig-caption="false"
  >
    <!-- HEADING -->
    <template #heading>
      <h1 v-if="titleLevel === '1'" v-html="title" />
      <h2 v-if="titleLevel === '2'" v-html="title" />
      <h3 v-if="titleLevel === '3'" v-html="title" />
    </template>
    <template #aboveExplanation>
      <div>
        <div
          v-for="reference in references"
          class="reference-text"
          :key="reference.num"
        >
          <ol>
            <span v-html="reference.authors" /> (<span v-html="reference.year" />). <a
              :href="reference.link"
              target="_blank"
            ><span v-html="reference.title" :class="{ report: reference.report }"/></a>
            <span v-if="reference.data_release">: U.S. Geological Survey data release</span>.
            <span v-if="reference.publication_info"> {{ reference.publication_info }}. </span>
            <span v-if="reference.journal">
              <span v-html="reference.journal_name" class="journal-name"></span>
              <span v-if="reference.journal_issue">, {{ reference.journal_issue }}</span>.
            </span>
            <span v-if="reference.doi" v-html="reference.doi"></span>
            <span v-else v-html="reference.link"></span>
          </ol>
        </div>
      </div>
    </template>
  </VizSection>
</template>

<script setup>
  import VizSection from '@/components/VizSection.vue';

  // define props
  defineProps({
    title: {
      type: String,
    },
    titleLevel: {
      type: String,
    },
    references: {
      type: Object,
    },
  })
</script>

<style scoped lang="scss">
  .reference-text {
    text-indent: -4rem;
    padding: 0 0 1rem;
  }
  .reference-text ol {
    list-style-type: none;
  }
  .journal-name {
    font-style: italic;
  }
  .report {
    font-style: italic;
  }
</style>