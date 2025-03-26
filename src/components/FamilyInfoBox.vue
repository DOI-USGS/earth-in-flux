<template>
  <div class="info-box">
    <figure v-if="activeFamily">
      <img :src="activeFamily.image" alt="Family illustration" />
      <figcaption v-html="activeFamily.caption"></figcaption>
    </figure>
    <p>{{ activeFamily.text }}</p>
  </div>
</template>

  <script setup>
  import { computed } from "vue";

  const props = defineProps({
    activeFamily: {
      type: Object,
      required: true
    }
    });

  const isPlaceholder = computed(() =>
    props.activeFamily.text?.includes("Click on a fish family")
  );

  </script>
  
  <style scoped>
  .info-box {
    background-color: #f8f9fa;
    padding: 2rem;
    border-radius: 12px;
    margin-bottom: 1.5rem;
    display: flex;
    gap: 4rem;
    align-items: center;
    box-shadow: 0 2px 8px rgba(0, 0, 0, 0.2);
  }
  
  img {
  width: 300px;     
  height: auto;    
  object-fit: contain; 
  border-radius: 8px;
  }

  figcaption {
    margin-top: 0.5rem;
    font-size: 1.15rem;
    color: #666;
    line-height: 1.2;
    max-width: 350px;
    text-align: center;
  }
  
  p {
    margin: 0 ;
    font-size: 1.25rem;
    line-height: 1.4;
  }

  .placeholder {
  opacity: 0.8;
  font-style: italic;
  }

  </style>
  