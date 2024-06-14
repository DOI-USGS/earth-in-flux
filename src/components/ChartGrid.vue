<template>
    <div id="chartGrid" class="padded">
		<ChartCard @click.enter="showSubPage(item.project, item.pageRoute)" v-for="(item, index) in chartContent" :key="index"
			:id="item.pageRoute"
			:src="item.img_src"
			:alt="item.alt"
			:project="item.project"
            :view="view"
		/>
    </div>
</template>

<script setup>
	import { useRouter } from 'vue-router'

	import ChartCard from '@/components/ChartCard.vue'
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

	function showSubPage(project, pageRoute) {
		const projectRoute = project.replace(/\s+/g, '-').toLowerCase();
		router.push({ name: 'SubPage', params: { projectRoute, pageRoute } })
	}

</script>

<style>

#chartGrid{
	margin: 50px auto;
	display: flex;
	justify-content: center;
	gap: 50px;
	flex-wrap: wrap;
	gap: 50px;
	width: 100%;
	max-width: 60%;
}
</style>