<template>
    <div id="chartGrid" class="padded">
        <ChartCard @click.enter="showSubPage(item.project, item.vizRoute)" v-for="(item, index) in sortedChartContent" :key="index"
            :id="item.vizRoute"
            :src="getThumb(item.img_src)"
            :alt="item.alt"
            :description="item.description"
        />
    </div>
</template>

<script setup>
    import { useRouter } from 'vue-router';
    import { computed, onMounted } from 'vue';
    import { isMobile } from 'mobile-device-detect';

    import ChartCard from '@/components/ChartCard.vue';
    import ChartGridContent from '@/assets/content/ChartGrid.js';

    const props = defineProps({
        view: {
            type: String,
            default: ``
        },
    })
    
    // global variables
    const router = useRouter();
    const mobileView = isMobile;
    const chartContent = ChartGridContent.chartGridItems;

    // set up filtered chart data as computed property
    const filteredChartContent = computed(() => {
        return props.view == 'all' ? chartContent : chartContent.filter(d => d.project.replace(/\s+/g, '-').toLowerCase() === props.view).sort((a,b) => (a.chartOrder > b.chartOrder) ? 1 : ((b.chartOrder > a.chartOrder) ? -1 : 0))
    });    
    
    // computed property for randomized chart content - only randomized on landing view
    const sortedChartContent = computed(() => {
        return props.view == 'all' ? shuffle([...filteredChartContent.value]) : filteredChartContent.value; // clone array to avoid mutating original
    });

    // Declare behavior on mounted
    // functions called here
    onMounted(async () => {
        if (mobileView) {
            const cards = document.querySelectorAll('.chart');

            window.addEventListener('scroll', function() {
                // add event on scroll
                cards.forEach(element => {
                    //for each .chart
                    if (isInViewport(element)) {
                        //if in Viewport
                        element.classList.add("expanded");
                    } else {
                        element.classList.remove("expanded");
                    }
                });
            }, false);
        }

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

    function showSubPage(project, vizRoute) {
        const projectRoute = project.replace(/\s+/g, '-').toLowerCase();
        router.push({ name: 'SubPage', params: { projectRoute, vizRoute } })
    }

    function getThumb(pic) {
        return 'https://labs.waterdata.usgs.gov/visualizations/thumbnails/'+pic
    }

    function isInViewport(elem) {
        const distance = elem.getBoundingClientRect();
        return (
            distance.top >= 0 &&
            distance.left >= 0 &&
            distance.bottom <= (window.innerHeight || document.documentElement.clientHeight) &&
            distance.right <= (window.innerWidth || document.documentElement.clientWidth)
        );
    };

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
        max-width: 100%;
    }
}
</style>