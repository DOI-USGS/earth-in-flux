<template>
  <VizSection
    id="authors"
    :figures="false"
    :fig-caption="false"
  >
    <!-- HEADING -->
    <template #heading>
      <hr>
      <h1 v-if="titleLevel === '1'" v-html="title" />
      <h2 v-if="titleLevel === '2'" v-html="title" />
      <h3 v-if="titleLevel === '3'" v-html="title" />
    </template>
    <template #aboveExplanation>
      <div
        v-if="showAuthors"
        id="author-container"
        class="text-content"
      >
        <p>
          This visualization was developed by the <a href='https://labs.waterdata.usgs.gov/visualizations/' target='_blank'>USGS Vizlab</a>
          <span id="primary-author-statment">
            and led by 
            <span
              v-for="(author, index) in authors" 
              :id="`initial-${author.initials}`"
              :key="`${author.initials}-attribution`"
              :class="'author first'"
            >
              <a
                :href="author.profile_link"
                target="_blank"
                v-text="author.fullName"
              />
              <span v-if="index != Object.keys(authors).length - 1 && Object.keys(authors).length > 2">, </span>
              <span v-if="index == Object.keys(authors).length - 2"> and </span>
            </span>.
          </span>
        </p>
      </div>
    </template>
  </VizSection>
</template>

<script setup>
  import { ref, onMounted } from 'vue';
  import VizSection from '@/components/VizSection.vue';

  // define props
  const props = defineProps({
    title: {
      type: String,
    },
    titleLevel: {
      type: String,
    },
    authors: {
      type: Object,
    },
  })

  // Turn on or off attribution for all authors
  const showAuthors = ref(null);

  onMounted(() => {
    showAuthors.value = props.authors.length > 0;
  });

</script>

<style>
  #authors {
    font-style: italic;
    font-weight: 300;
  }
</style>