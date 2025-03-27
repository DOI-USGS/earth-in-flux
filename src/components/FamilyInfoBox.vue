<template>
  <div class="info-box" v-if="activeFamily">
    
      <!-- FAMILY LEVEL -->
    <div v-if="activeFamily.type === 'family' || !activeFamily.type" class="family-content">
      <figure>
        <img :src="activeFamily.image" alt="Fish silhouette" />
        <figcaption v-if="activeFamily.caption" v-html="activeFamily.caption"></figcaption>
      </figure>
      <div class="family-text">
        <h2 v-if="activeFamily.name">{{ activeFamily.name }}</h2>
        <p v-if="activeFamily.economicValue"><strong>{{ activeFamily.economicValue }}</strong> in total economic value</p>
        <p v-if="activeFamily.speciesCount"><strong>{{ activeFamily.speciesCount }}</strong> species that are recreationally fished</p>
        <p>{{ activeFamily.text }}</p>
      </div>
    </div>
    <!-- SPECIES LEVEL -->
    <div v-else-if="activeFamily.type === 'species'" class="species-content">
      <div class="species-header">
        <img :src="activeFamily.image" alt="Fish silhouette small" class="species-icon" />
        <span class="family-label">{{ activeFamily.family }}</span>
      </div>
      <hr />
      <h2>{{ activeFamily.name }}</h2>
      <p class="species-subtitle">Species</p>
      <p><strong>{{ activeFamily.economicValue }}</strong> in total economic value</p>
      <p><strong>{{ activeFamily.countryCount }}</strong> countries where it is recreationally fished</p>
    </div>

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

  </script>
  
  <style scoped>
 .info-box {
  background-color: #f8f9fa;
  padding: 2rem;
  border-radius: 12px;
  margin-bottom: 2rem;
  margin-top: 2rem;
  box-shadow: 0 2px 8px rgba(0, 0, 0, 0.2);
  display: block;
}

.family-content {
  display: flex;
  align-items: center;
  gap: 2rem;
}

.family-text {
  max-width: 600px;
   font-size: 1.4rem;
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
  text-align: center;
}

.species-content {
  display: flex;
  flex-direction: column;
  font-size: 1.4rem
}

.species-header {
  display: flex;
  align-items: center;
  gap: 1rem;
}

.species-icon {
  width: 100px;
}

.family-label {
  font-size: 1.4rem;
  color: #333;
}

.species-subtitle {
  font-size: 1.4rem;
  color: #666;
  margin-bottom: 0.5rem;
}

.species-content hr {
  margin-top: 1rem;
  margin-bottom: 0.5rem
}


  </style>
  