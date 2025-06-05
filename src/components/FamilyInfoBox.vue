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
        <p v-if="activeFamily.type === 'family'" class="species-subtitle">Family</p>
        <p v-if="activeFamily.economicValue">
          <strong>{{ activeFamily.economicValue }}</strong> total consumptive use value
        </p>
        <p v-if="activeFamily.speciesCount">
          <strong>{{ activeFamily.speciesCount }}</strong>
          {{
            activeFamily.speciesCount === 1
              ? 'species that is recreationally harvested for human consumption'
              : 'species that are recreationally harvested for human consumption'
          }}
        </p>
        <p>{{ activeFamily.text }}</p>
      </div>
    </div>
    <!-- OTHER-SPECIES LUMPING BIN -->
    <div v-else-if="activeFamily.name === 'Other'" class="species-content">
      <div class="species-header">
        <img :src="activeFamily.image" alt="Fish silhouette small" class="species-icon" />
        <span class="family-label">{{ activeFamily.family }}</span>
      </div>
      <hr />
      <h2>Other</h2>
      <p class="species-subtitle">
        <strong>{{ activeFamily.lumpedSpeciesCount }}</strong>
        {{ activeFamily.lumpedSpeciesCount === 1 ? 'species' : 'species' }}
      </p>
      <p>
        <strong>${{ roundValue(activeFamily.economicValue) }}</strong> total consumptive use value
      </p>
      <p>
        <strong>{{ activeFamily.countryCount }}</strong>
        {{
          activeFamily.countryCount === 1
            ? 'country where they are recreationally harvested for human consumption'
            : 'countries where they are recreationally harvested for human consumption'
        }}
      </p>
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
      <p>
        <strong>{{ activeFamily.economicValue }}</strong> total consumptive use value
      </p>
      <p>
        <strong>{{ activeFamily.countryCount }}</strong>
        {{
          activeFamily.countryCount === 1
            ? 'country where it is recreationally harvested for human consumption'
            : 'countries where it is recreationally harvested for human consumption'
        }}
      </p>
    </div>
  </div>
</template>

<script setup>
defineProps({
  activeFamily: {
    type: Object,
    required: true
  }
})
function roundValue(value) {
  const number = Number(value)
  if (isNaN(number)) return '0'
  return Math.round(number).toLocaleString()
}
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
  font-size: 1.4rem;
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
  font-weight: 700;
}

.species-subtitle {
  font-size: 1.4rem;
  color: #666;
  margin-bottom: 0.5rem;
  margin-top: -1rem;
}

.species-content hr {
  margin-top: 1rem;
  margin-bottom: 0.5rem;
}

@media (max-width: 700px) {
  .info-box {
    padding: 1rem;
  }

  .family-content figure {
    align-items: center;
    justify-content: center;
  }

  .family-text {
    max-width: 100%;
    font-size: 1.2rem;
  }

  img {
    width: 100%;
    max-width: 400px;
    height: auto;
    align-self: center;
  }

  figcaption {
    text-align: center;
  }

  .species-content {
    font-size: 1.2rem;
  }

  .species-icon {
    width: 60px;
  }

  .family-label {
    font-size: 1.2rem;
  }

  .species-subtitle {
    font-size: 1.2rem;
  }
}
</style>
