var awsIot = require('aws-iot-device-sdk');


var device = awsIot.device({
   keyPath: __dirname+'/testdago.private.key',
  certPath: __dirname+'/testdago.cert.pem',
    caPath: __dirname+'/root-CA.crt',
  clientId: 'Nanopore_2',
      host: 'a1awxf1oj7f69s.iot.us-east-2.amazonaws.com'
});

device
  .on('connect', function() {

device.subscribe('air_filter/nanopore_2');
