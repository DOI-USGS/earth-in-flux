<template>
  <div class="info-box" v-if="activeCountry">
    <div class="country-content">
      <div class="country-text">
        <div
          v-if="!activeCountry.name"
          id="intro-content-container"
        >
          <img class="map-image" :src="getImageURL('fish_as_food_continent_map.png')" alt="">
          <p>{{ activeCountry.text }}</p>
        </div>
        <div v-if="activeCountry.name">
          <h2>{{ activeCountry.name }}</h2>
          <hr>
          <span
            class="card-line"
          >
            <img class="icon-image" :src="getSvgURL('noun-fishermen-305219')" alt="">
            <span><b>{{ activeCountry.n_fishers }}</b> recreational fishers ({{ activeCountry.participation_rate }} % of population)</span>
          </span>
          <span
            class="card-line"
          >
            <img class="icon-image" :src="getSvgURL('noun-net-pen-1559229')" alt="">
            <span><b>{{ activeCountry.consum_harv_kg }} kg</b> of consumable harvest</span>
          </span>
          <span>
            <ul
              id="guild-breakdown-container"
            >
              <li
                v-for="guild, index in guildSummary"
                :key="index"
              >
                <span
                  class="card-line"
                >
                  <img class="icon-image" :src="identifySvgURL(guild.guild)" alt="">
                  <span v-if="guild.guild != 'unknown'"><b>{{ guild.percent }}%</b> {{ guild.guild }}-water species</span>
                  <span v-if="guild.guild == 'unknown'"><b>{{ guild.percent }}%</b> species with unknown thermal guild</span>
                </span>
              </li>
            </ul>
          </span>
          <span
            class="card-line"
          >
            <img class="icon-image" :src="getSvgURL('noun-global-warming-7927954')" alt="">
            <span><b>{{ activeCountry.mean_vul }} out of 1</b> average <a href='https://onlinelibrary.wiley.com/doi/10.1111/gcb.15768?af=R' target='_blank'>climate vulnerability score</a> for harvested species</span>
          </span>
        </div>
      </div>
    </div>
  </div>
</template>

<script setup>
import { computed } from "vue";

const props = defineProps({
  activeCountry: {
    type: Object,
    required: true
  }
})

const guildSummary = computed(() => {
  console.log(props.activeCountry.harv_breakdown)
  return props.activeCountry.harv_breakdown.filter(d => d.percent)
});

function getSvgURL(filename) {
  return new URL(`../assets/svgs/${filename}.svg`, import.meta.url).href
}

function getImageURL(filename) {
  return new URL(`../assets/images/${filename}`, import.meta.url).href
}

function identifySvgURL(guild) {
  let filename = '';
  switch (guild) {
    case 'warm':
      console.log('it is warm')
      filename = 'noun-water-temperature-743356'
      return new URL(`../assets/svgs/${filename}.svg`, import.meta.url).href
      // break;
    case 'cool':
      filename = 'noun-water-temperature-743353';
      return new URL(`../assets/svgs/${filename}.svg`, import.meta.url).href
    case 'cold':
      filename = 'noun-water-temperature-743352';
      return new URL(`../assets/svgs/${filename}.svg`, import.meta.url).href
    default:
      filename = 'noun-question-mark-7853683';
      return new URL(`../assets/svgs/${filename}.svg`, import.meta.url).href
  }

}

</script>

<style scoped>
.info-box {
    max-width: 60vw;
    margin-left: auto;
    margin-right: auto;
    @media screen and (min-width: 2000px) {
      max-width: 1000px;
    }
    @media screen and (max-width: 600px) {
        width: 100%;
        height: auto;
        min-width: auto;
        max-width: 100%;
    }
}
#intro-content-container {
  display: flex;
  flex-direction: row;
  align-items: center;
  gap: 2rem;
  @media screen and (max-width: 600px) {
    flex-direction: column;
  }
}
.map-image {
  width: 300px;
}
.country-content {
    /* min-width: 60vw; */
    display: flex;
    flex-direction: row;
    @media screen and (max-width: 600px) {
        flex-direction: column;
    }
}
.country-content h2 {
    padding-top: 0rem;
}
.country-content p {
    padding-bottom: 0.2rem;
}
.country-content hr {
  margin-top: 0.1rem;
  margin-bottom: 0.5rem;
}
.card-line {
  display: flex;
  flex-direction: row;
  flex-wrap: nowrap;
  align-items: center;
  margin-bottom: 1rem;
}
.icon-image {
  max-width: 50px;
  margin-right: 10px;
}
#guild-breakdown-container {
  list-style-type: none;
  display: flex;
  flex-direction: row;
  align-items: center;
  @media screen and (max-width: 600px) {
    flex-direction: column;
    align-items: start;
  }
}
#guild-breakdown-container li {
  margin-right: 1rem;
  @media screen and (max-width: 600px) {
    margin-right: 0rem;
  }
}
#guild-breakdown-container li:first-of-type {
  padding-top: 0rem;
}
#guild-breakdown-container li {
  padding-bottom: 0rem;
}
</style>
