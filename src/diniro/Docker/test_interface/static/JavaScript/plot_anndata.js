var ctx = document.getElementById('myChart').getContext('2d');
var scatterChart = new Chart(ctx, {
    type: 'scatter',
    data: {
        datasets: [{
            label: 'Scatter Dataset',
            data: pp,
            pointBackgroundColor: cc,
            pointBorderWidth : 0.8,
            pointBorderColor : 'black',
            pointRadius : '5'

        }]
    },
    options: {
        scales: {
            xAxes: [{
                type: 'linear',
                position: 'bottom'
            }]
        }
    }
});