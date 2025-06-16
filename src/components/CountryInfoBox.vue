<template>
  <div class="info-box" v-if="activeCountry">
    <div class="country-content">
      <div class="country-text">
        <div v-if="!activeCountry.name">
          <p>{{ activeCountry.text }}</p>
        </div>
        <div v-if="activeCountry.name">
          <h2>{{ activeCountry.name }}</h2>
          <hr>
          <p><b>Population:</b> <span>{{ activeCountry.population }}</span></p>
          <p><span>{{ activeCountry.n_fishers }}</span> recreational fishers (<span>{{ activeCountry.participation_rate }}</span>% of population)</p>
          <p><b>Total consumable harvest:</b> <span>{{ activeCountry.consum_harv_kg }}</span> kg</p>
          <p>Average <a href='https://onlinelibrary.wiley.com/doi/10.1111/gcb.15768?af=R' target='_blank'>climate vulnerability score</a> (0-1) of inland species harvested for consumption: <span>{{ activeCountry.mean_vul }}</span></p>
          <p
           v-for="guild, index in guildSummary"
           :key="index"
          >
            {{ guild.guild }} : {{ guild.percent }} % of harvest
          </p>
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
  return props.activeCountry.harv_breakdown.filter(d => d.percent)
});

</script>

<style scoped>
.info-box {
    width: 60vw;
    margin-left: auto;
    margin-right: auto;
    height: 200px;
    @media screen and (max-width: 600px) {
        width: 100%;
        height: 400px;
    }
}
.country-content {
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
</style>
