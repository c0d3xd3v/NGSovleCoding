import QtQuick
import QtQuick.Controls

ToolButton {
    id: toolButton
    text: qsTr("back")
    onClicked: function() {
        console.log(swipeView.currentIndex)
        swipeView.currentIndex = 0
    }
}
