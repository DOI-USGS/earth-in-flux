<template>
    <div id="chartGrid">
        <ChartIcon 
            :id="projectRoute"
            :vizRoute=null
            :projectRoute="projectRoute"
        />
        <ChartIcon v-for="item in filteredChartContent" :key="item.vizRoute"
            :id="item.vizRoute"
            :src="getThumb(item.img_src)"
            :alt="item.alt"
            :active="item.vizRoute === vizRoute"
            :vizRoute="item.vizRoute"
            :projectRoute="projectRoute"
        />
    </div>
</template>

<script setup>
    import { computed } from 'vue';

    import ChartIcon from '@/components/ChartIcon.vue';
    import ChartGridContent from '@/assets/content/ChartGrid.js';

    const props = defineProps({
        projectRoute: {
            type: String,
            default: ``
        },
        vizRoute: {
            type: String,
            default: ``
        },
    })
    
    // global variables
    const chartContent = ChartGridContent.chartGridItems;

    // set up filtered chart data as computed property
    const filteredChartContent = computed(() => {
        return chartContent.filter(d => d.project.replace(/\s+/g, '-').toLowerCase() === props.projectRoute).sort((a,b) => (a.chartOrder > b.chartOrder) ? 1 : ((b.chartOrder > a.chartOrder) ? -1 : 0))
    });    

    function getThumb(pic) {
        return 'https://labs.waterdata.usgs.gov/visualizations/thumbnails/'+pic
    }

</script>

<style>
#chartGrid{
    margin: 8rem auto 2rem auto;
    display: flex;
    justify-content: center;
    gap: 50px;
    flex-wrap: wrap;
    width: 100%;
    max-width: 70rem;
    @media screen and (max-width: 600px) {
        gap: 25px;
        max-width: 95%;
    }
}
</style>