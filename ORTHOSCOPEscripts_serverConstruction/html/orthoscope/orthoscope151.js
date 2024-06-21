//function checkradio( disp ) {
//   document.getElementById('id_display1').style.display = disp;
//   document.getElementById('id_display2').style.display = disp;
//   document.getElementById('id_display3').style.display = disp;
//   document.getElementById('id_display4').style.display = disp;
//   document.getElementById('id_display5').style.display = disp;
//}

function checkradioSeqChr( disp ) {
   document.getElementById('id_displaySeqChr1').style.display = disp;
   document.getElementById('id_displaySeqChr2').style.display = disp;
   document.getElementById('id_displaySeqChr3').style.display = disp;
}

$(function(){
  $('.checkAll').on('change', function() {
    $('input[name="taxon_checkbox"]').prop('checked', this.checked);
    $('input[name="taxon_checkbox_focal"]').prop('checked', this.checked);
  });
});

$(function(){
  $('.checkAll_focal').on('change', function() {
    $('input[name="taxon_checkbox_focal"]').prop('checked', this.checked);
  });
});


//$(function() {
//  // ボタンをクリックしたら発動
//  $('#btn').click(function() {
//    if ($('input[name="taxon_checkbox"]').prop('checked') || $('input[name="taxon_checkbox_focal"]').prop('checked')) {
//      $('input[name="taxon_checkbox"]').prop('checked', false);
//      $('input[name="taxon_checkbox_focal"]').prop('checked', false);
//    } else {
//      $('input[name="taxon_checkbox"]').prop('checked', true);
//      $('input[name="taxon_checkbox_focal"]').prop('checked', true);
//    }
//  });
//});


//$(function() {
//  // ボタンをクリックしたら発動
//  $('#btn').click(function() {
//    if ($('input[name="taxon_checkbox"]').prop('checked') || $('input[name="taxon_checkbox_focal"]').prop('checked')) {
//      $('input[name="taxon_checkbox"]').prop('checked', false);
//      $('input[name="taxon_checkbox_focal"]').prop('checked', false);
//    } else {
//      $('input[name="taxon_checkbox"]').prop('checked', true);
//      $('input[name="taxon_checkbox_focal"]').prop('checked', true);
//    }
//  });
//});


$(function() {
  $('#btn_focal').click(function() {
    if ($('input[name="taxon_checkbox_focal"]').prop('checked')) {
      $('input[name="taxon_checkbox_focal"]').prop('checked', false);
    } else {
      $('input[name="taxon_checkbox_focal"]').prop('checked', true);
    }
  });
}); 


$(function() {
  $('#my_form').submit('click', function(e) {
    e.preventDefault();

    //var consoleInfo = console.info($('#my_form').get(0).submit);
    //console.log(consoleInfo)

    var fd = new FormData($('#my_form').get(0));

    //if (fd.get("input_file").size == 0) {
    //    fd.delete("input_file"); 
    //}

    if (fd.get("input_treeFile").size == 0) {
        fd.delete("input_treeFile");
    }

    //for (var key of fd.keys()) {
    //    console.log(key);
    // }

/*
    for (item of fd){
      console.log(item)
    }
*/

    $(document)
    .ajaxStart(function(){
      //$('#prog').show();
      $('#result').html('<img src="orthoscope.gif" alt="" />');
    })
    console.log("test1")
    $.ajax({
      url: '/cgi-bin/orthoscope151.py',
      type: 'post',
      dataType: 'text',
      data: fd,
      processData: false,
      contentType: false,
    })
    .done(function(response) {
      console.log("testDone")
      $('#result').html(response);
    })
    $(document)
    .ajaxError(function(e, xhr, opts, err){
      console.log("testError_timeout")
      $('#result').html('Error: timeout' + err);
    })
    //.fail(function() {
    //  $('#result').html('Failed.');
    //})
    //$(document)
    //.ajaxStop(function(){
    //  $('#prog').hide();
    //})
  });
  return false;
});
