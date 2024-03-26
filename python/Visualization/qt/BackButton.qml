import QtQuick
import QtQuick.Controls

ToolButton {
    id: toolButton
    text: qsTr("back")
    onClicked: function() {
        swipeView.currentIndex = 0
    }
}
