$(document).ready(function(){
  $("button").click(function(){
    $.ajax({url: 'about/',
        type : 'GET',
        data: { test: 'ajaxmethod'},
        success: function(data){

    }});
  });
});