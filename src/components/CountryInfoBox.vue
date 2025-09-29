<template>
  <div class="info-box" v-if="activeCountry">
    <div class="country-content">
      <div class="country-text">
        <div
          v-if="!activeCountry.name"
          id="intro-content-container"
        >
          <img class="map-image" :src="getImageURL('fish_as_food_continent_map.png')" alt="">
          <span>
            <p>{{ activeCountry.text }}</p>
            <br>
            <p
              v-if="mobileView"
            >
              {{ activeCountry.prompt_mobile }}
            </p>
            <p
              v-else
            >
              {{ activeCountry.prompt_desktop }}
            </p>
          </span>
        </div>
        <div v-if="activeCountry.name"
          class="country-box"
          :class="activeCountry.continent"
        >
          <div
            id="country-map-container"
          >
            <h2>{{ activeCountry.name }}</h2>
            <img class="map-image-small" :src="getImageURL('fish_as_food_continent_map.png')" alt=""></img>
          </div>
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
            <span><b>{{ activeCountry.consum_harv_kg }} kg</b> total consumable harvest<sup>1</sup></span>
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
                  <img class="icon-image guild" :src="identifySvgURL(guild.guild)" alt="">
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
import { isMobile } from 'mobile-device-detect';

const props = defineProps({
  activeCountry: {
    type: Object,
    required: true
  }
})

// global variables
const mobileView = isMobile;

const guildSummary = computed(() => {
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
    width: 100%;
    margin-left: auto;
    margin-right: auto;
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
  min-height: 285px;
  @media screen and (max-width: 1000px) {
    flex-direction: column;
  }
}
.map-image {
  width: 300px;
}
.map-image-small {
  width: 100px;
}
.country-content h2 {
  padding: 0;
  line-height: 1;
}
.country-content p {
  padding-bottom: 0.2rem;
}
.country-content hr {
  margin-top: 0.1rem;
  margin-bottom: 0.5rem;
}
.country-box {
  border-left: 12px black solid;
  padding-left: 10px;
  border-radius: 8px;
}
.Europe {
  border-color: var(--mid-brown);
}
.Africa {
  border-color: var(--mid-teal);
}
.Oceania {
  border-color: var(--mid-purple);
}
.North_America {
  border-color: var(--mid-red);
}
.Asia {
  border-color: var(--dark-blue);
}
.South_America {
  border-color: var(--grey-blue);
}
#country-map-container {
  display: flex;
  flex-direction: row;
  justify-content: space-between;
  align-items: center;
  padding: 0.2rem 0 0.2rem 0;
}
.card-line {
  display: flex;
  flex-direction: row;
  flex-wrap: nowrap;
  align-items: center;
  margin-bottom: 1rem;
}
.card-line:last-of-type {
  margin-bottom: 0rem;
}
.icon-image {
  max-width: 50px;
  margin-right: 10px;
}
.icon-image.guild {
  max-width: 45px;
  margin-right: 10px;
}
#guild-breakdown-container {
  list-style-type: none;
  display: flex;
  flex-direction: row;
  align-items: center;
  max-width: 100%;
  font-size: 1.6rem;
  @media screen and (max-width: 600px) {
    flex-direction: column;
    align-items: start;
    max-width: 100%;
  }
}
ul#guild-breakdown-container {
  padding-inline-start: 60px;
  @media screen and (max-width: 600px) {
    padding-inline-start: 40px;
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
