<template>
  <VizSection
    id="references"
    :figures="false"
    :fig-caption="false"
  >
    <!-- HEADING -->
    <template #heading>
      <h1
        v-if="titleLevel === '1'"
        v-html="title"
      />
      <h2
        v-if="titleLevel === '2'"
        v-html="title"
      />
      <h3
        v-if="titleLevel === '3'"
        v-html="title"
      />
    </template>
    <template #aboveExplanation>
      <div>
        <div
          v-for="reference in references"
          :key="reference.num"
          class="reference-text"
        >
          <ol>
            <span v-html="reference.authors" /> (<span v-html="reference.year" />). <a
              :href="reference.link"
              target="_blank"
            ><span
              :class="{ report: reference.report }"
              v-html="reference.title"
            /></a>
            <span v-if="reference.data_release">: U.S. Geological Survey data release</span>.
            <span v-if="reference.publication_info"> {{ reference.publication_info }}. </span>
            <span v-if="reference.journal">
              <span
                class="journal-name"
                v-html="reference.journal_name"
              />
              <span v-if="reference.journal_issue">, {{ reference.journal_issue }}</span>.
            </span>
            <span
              v-if="reference.doi"
              v-html="reference.doi"
            />
            <span
              v-else
              v-html="reference.link"
            />
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
      default: ""
    },
    titleLevel: {
      type: String,
      default: ""
    },
    references: {
      type: Object,
      default: () => ({})
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