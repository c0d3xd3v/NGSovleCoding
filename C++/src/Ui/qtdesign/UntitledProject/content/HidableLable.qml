import QtQuick 2.15
import QtQuick.Controls 2.15

Label {
    text: "label"
    visible: true
    function test() {
    }
    Component.onCompleted: test()
}
