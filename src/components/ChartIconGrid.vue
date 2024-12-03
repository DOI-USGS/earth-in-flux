<template>
    <div id="chartGrid" class="padded">
        <ChartIcon v-for="item in filteredChartContent" :key="item.vizRoute"
            :id="item.vizRoute"
            :src="getThumb(item.img_src)"
            :alt="item.alt"
            :description="item.description"
            :active="item.vizRoute === viz"
            :vizRoute="item.vizRoute"
        />
    </div>
</template>

<script setup>
    import { computed } from 'vue';

    import ChartIcon from '@/components/ChartIcon.vue';
    import ChartGridContent from '@/assets/content/ChartGrid.js';

    const props = defineProps({
        project: {
            type: String,
            default: ``
        },
        viz: {
            type: String,
            default: ``
        },
    })
    
    // global variables
    const chartContent = ChartGridContent.chartGridItems;

    // set up filtered chart data as computed property
    const filteredChartContent = computed(() => {
        return chartContent.filter(d => d.project.replace(/\s+/g, '-').toLowerCase() === props.project).sort((a,b) => (a.chartOrder > b.chartOrder) ? 1 : ((b.chartOrder > a.chartOrder) ? -1 : 0))
    });    

    function getThumb(pic) {
        return 'https://labs.waterdata.usgs.gov/visualizations/thumbnails/'+pic
    }

</script>

<style>
#chartGrid{
    margin: 50px auto;
    display: flex;
    justify-content: center;
    gap: 50px;
    flex-wrap: wrap;
    width: 100%;
    max-width: 60%;
    @media screen and (max-width: 600px) {
        gap: 25px;
        max-width: 95%;
    }
}
</style>