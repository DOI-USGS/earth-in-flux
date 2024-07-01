<template>
    <div id="chartGrid" class="padded">
        <ChartCard @click.enter="showSubPage(item.project, item.vizRoute)" v-for="(item, index) in randomizedChartContent" :key="index"
            :id="item.vizRoute"
            :src="getThumb(item.img_src)"
            :alt="item.alt"
        />
    </div>
</template>

<script setup>
    import { useRouter } from 'vue-router';
    import { computed } from 'vue';

    import ChartCard from '@/components/ChartCard.vue';
    import ChartGrid from '@/assets/content/ChartGrid.js';

    const props = defineProps({
        view: {
            type: String,
            default: ``
        },
    })
    
    // global variables
    const router = useRouter()
    const chartContent = ChartGrid.chartGridItems;

    // set up filtered chart data as computed property
    const filteredChartContent = computed(() => {
        return props.view == 'all' ? chartContent : chartContent.filter(d => d.project.replace(/\s+/g, '-').toLowerCase() === props.view)
    });

    // function to shuffle an array
    function shuffle(array) {
        let currentIndex = array.length,  randomIndex;

        // While there are more elements to shuffle
        while (currentIndex != 0) {

            // Picks a remaining element
            randomIndex = Math.floor(Math.random() * currentIndex);
            currentIndex--;

            // Swaps it with the current element
            [array[currentIndex], array[randomIndex]] = [
                array[randomIndex], array[currentIndex]];
        }

        return array;
    }

    // computed property for randomized chart content
    const randomizedChartContent = computed(() => {
        return shuffle([...filteredChartContent.value]); // clone array to avoid mutating original
    });

    function showSubPage(project, vizRoute) {
        const projectRoute = project.replace(/\s+/g, '-').toLowerCase();
        router.push({ name: 'SubPage', params: { projectRoute, vizRoute } })
    }

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
}
</style>